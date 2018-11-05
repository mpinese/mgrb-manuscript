import math
import memfiles
import parseopt
import parseutils
import random
import sets
import sequtils
import strutils
import tables
import terminal

const VERSION = "0.0.4, 13 Mar 2018"
const AUTHOR = "Mark Pinese <m.pinese@garvan.org.au>"


type
  Dosages = tuple[
    samples: seq[string],
    vids: seq[string],
    vid2idx: Table[string, int],
    afs: seq[float],
    afbins: seq[int],
    afbin2idx: seq[seq[int]],
    dosages: seq[int8]
  ]
  VarCoef = tuple[af: float, afbin: int, coef: float]
  Model = tuple[
    id: string,
    offset: float,
    coefs: Table[string, VarCoef]   # Keyed by vid
  ]
  Models = Table[string, Model]     # Keyed by id
  SummaryStatistics = array[4, float]
  SamplingReference {.pure.} = enum
    Internal, External
  OutputFormat {.pure.} = enum
    Complete, Summary


proc scanDosages(f: MemFile): tuple[samples: seq[string], nvariants:int] = 
  # Count lines and columns in the file.  Check for sample uniqueness.
  var
    header_read = false
    nvariants = 0
    samples = @[""]
  for slice in f.memSlices:
    if header_read == false:
      header_read = true
      samples = ($slice).strip(leading=false).split(sep="\t")[1..^1]
      var sample_set = initSet[string](sets.rightSize(samples.len))
      for sample in result.samples:
        if sample_set.containsOrIncl(sample) == true:
          raise newException(Exception, "Sample " & sample & " present more than once in dosages file")
    else:
      nvariants += 1
  result = (samples, nvariants)


proc loadDosages(path: string, n_afbins=50): Dosages = 
  stderr.write("Loading dosages from " & path & ", " & $n_afbins & " AF bins...\n")

  # Scan the file to get sample IDs and variant count.
  # Verify sample ID uniqueness.
  let
    f = memfiles.open(path, mode=fmRead)
    (samples, nvariants) = scanDosages(f)
    nsamples = samples.len
  result.samples = samples
  stderr.write("  " & $nvariants & " variants x " & $nsamples & " samples found.  Allocating...")

  # Preallocate data
  result.dosages = newSeq[int8](nsamples*nvariants)
  result.vids = newSeq[string](nvariants)
  result.afs = newSeq[float](nvariants)
  result.afbins = newSeq[int](nvariants)
  result.vid2idx = initTable[string, int](initialSize=tables.rightSize(nvariants))
  result.afbin2idx = newSeqWith(n_afbins, newSeq[int]())
  stderr.write("  Reading...\n")

  # Read the data and verify vid uniqueness.
  var 
    i = 0
    vid_set = initSet[string](sets.rightSize(nvariants))
    header_read = false
    buffer: TaintedString = ""
    vid: string
    dosage: char
    k: int
    nmissing: int
    nalt: int
    ndiscarded = 0
  
  for line in f.lines(buffer):
    if header_read == false:
      header_read = true
      continue

    if i %% 10000 == 0:
      stderr.eraseLine()
      stderr.write("    " & $i & " / " & $nvariants & " variants")
    
    k = parseUntil(line, vid, {'\t'}, 0) + 1

    if line.len - k < nsamples:
      raise newException(Exception, "Insufficient genotypes for variant " & vid & ".  Expected " & $nsamples & ", found at most " & $(line.len - k))

    if vid_set.containsOrIncl(vid) == true:
      raise newException(Exception, "Variant " & vid & " present more than once in dosages file")

    let dosages_offset = i*nsamples
    nmissing = 0
    nalt = 0
    for j in 0..<nsamples:
      dosage = line[k+j]
      if not (dosage == '.' or dosage == '0' or dosage == '1' or dosage == '2'):
        raise newException(Exception, "Invalid dosage value " & $dosage & " found in dosage file, variant " & vid & " sample " & result.samples[j])
      if dosage == '.':
        result.dosages[dosages_offset + j] = -1
        nmissing += 1
      else:
        let dosage_int = int8(ord(dosage) - ord('0'))
        result.dosages[dosages_offset + j] = dosage_int
        nalt += dosage_int

    if nmissing.float / nsamples.float > 0.01:
      # stderr.write("    Too many genotypes missing for variant " & vid & ", discarding\n")
      # Simulate a variant discard by setting vid to "" (disallowed in loadModels), and
      # not adding this variant to the AF bin array.  The variant will then not be 
      # used for score calculation (vid == ""), or be included in resampled null models
      # (not in AF bins).
      vid = ""
      result.afs[i] = NaN
      result.afbins[i] = -1
      ndiscarded += 1
    else:
      # Add the variant to all data structures.
      result.vid2idx[vid] = i
      result.afs[i] = nalt.float / (2*(nsamples - nmissing)).float
      result.afbins[i] = min(n_afbins - 1, floor(result.afs[i] * n_afbins.float).int)
      result.afbin2idx[result.afbins[i]].add(i)
    result.vids[i] = vid
    i += 1

  # Check AF bin occupancy
  var min_occupancy = result.afbin2idx[0].len
  for i in 1..<n_afbins:
    min_occupancy = min(min_occupancy, result.afbin2idx[i].len)

  stderr.write("\n  Loaded " & $(nvariants - ndiscarded) & " variants x " & $nsamples & " samples, smallest AF bin size: " & $min_occupancy & "; " & $ndiscarded & " variants discarded.\n")


proc loadModels(path: string, n_afbins: int, dosages: Dosages): Models = 
  stderr.write("Loading models from " & path & ", " & $n_afbins & " AF bins...")
  result = initTable[string, Model]()
  let f = system.open(path, mode=fmRead)
  discard f.readLine()
  for line in f.lines:
    let
      fields = line.strip(leading=false).split(sep="\t")
      model_id = fields[0]
      vid = fields[1]
      coef = if fields[3] == "NA": 0.0 else: fields[3].parseFloat

    if vid == "OFFSET":
      result[model_id].offset = coef
      continue

    if vid == "":
      raise newException(Exception, "Invalid missing variant id encountered.")

    let
      af = fields[2].parseFloat
      afbin = min(n_afbins - 1, floor(af * n_afbins.float).int)
      dosages_af = if dosages.vid2idx.hasKey(vid): dosages.afs[dosages.vid2idx[vid]] else: af
    if af - dosages_af > 0.05 or af - dosages_af < -0.05:
      stderr.write("\n  " & model_id & " AF mismatch for vid " & vid & ": Dosages " & $dosages_af & ", Model " & $af)
    if not result.hasKey(model_id):
      result[model_id] = (id:model_id, offset:0.0, coefs:initTable[string, VarCoef]())

    result[model_id].coefs[vid] = (af:af, afbin:afbin, coef:coef)

  stderr.write("\n  Loaded " & $result.len & " models\n")


proc calcValues(model: Model, dosages: Dosages): seq[float] = 
  let n_samples = dosages.samples.len
  result = newSeq[float](n_samples)

  for i in 0..<n_samples:
    result[i] = model.offset

  for vid, coef in model.coefs.pairs:
    if dosages.vid2idx.hasKey(vid):
      let dosages_offset = dosages.vid2idx[vid]*n_samples
      for i in 0..<n_samples:
        if dosages.dosages[dosages_offset + i] == -1:
          # Missing genotype for this variant:sample.
          # Impute with empirical allele dosage based on samples with non-missing genotypes.
          # Variant allele frequency is dosages.afs[dosages.vid2idx[vid]]
          result[i] += 2.0*dosages.afs[dosages.vid2idx[vid]]*coef.coef
        else:
          # Genotype available for this variant:sample.
          # Variant allele dosage is dosages.dosages[dosages_offset + i]
          result[i] += float(dosages.dosages[dosages_offset + i])*coef.coef
    else:
      # Variant is missing from the dosages.  Impute with the model-reported
      # allele frequency (can't use empirical frequency from dosages as the
      # variant's not present).  The effect will be an equal translation of
      # the score for all individuals, so relative differences will be 
      # preserved.
      for i in 0..<n_samples:
        result[i] += 2.0*coef.af*coef.coef


proc generateNullModel(model: Model, dosages: Dosages, af_source: SamplingReference, seed: int): Model = 
  randomize(seed)

  result.id = model.id
  result.offset = model.offset
  result.coefs = initTable[string, VarCoef](initialSize=tables.rightSize(model.coefs.len))

  for vid, coef in model.coefs.pairs:
    if dosages.vid2idx.hasKey(vid):
      # Select a new variant from dosages with matching allele frequency.
      # There are two choices for target AF:
      #   - empirical AF in the dosages (SamplingReference.Internal)
      #   - AF as reported in the model (SamplingReference.External)
      # 
      # The first case behaves similarly to a permutation test, and 
      # effectively tests for a difference in score distribution shape between
      # the true scores and the resampled.  Essentially, are the true
      # scores more/less dispersed than expected by chance?  Mean shifts
      # will not be tested.
      #
      # The second case tests for mean shift caused by AF differences
      # between the dosages file and the AFs reported in the model file,
      # as well as potentially distribution shape differences.
      let afbin = 
        if af_source == SamplingReference.Internal:
          dosages.afbins[dosages.vid2idx[vid]]
        else:  # af_source = SamplingReference.External
          coef.afbin

      # Note: sampling with replacement.  Shouldn't matter in almost all cases.
      # TODO: Ideally should not select variants in LD.  Difficult to implement though.
      # One rough approach could be to enforce a minimum distance.
      # This sampling with replacement *will* be an issue for WGP.
      let
        new_vid_idx = dosages.afbin2idx[afbin][random(dosages.afbin2idx[afbin].len)]
        new_vid = dosages.vids[new_vid_idx]
      result.coefs[new_vid] = (af:dosages.afs[new_vid_idx], afbin:afbin, coef:coef.coef)
    else:
      # This locus was absent from dosages, so its null equivalent should
      # be missing too.  Easily done by leaving it alone.
      result.coefs[vid] = coef


proc calcSummaryStatistics(values: seq[float]): SummaryStatistics = 
  result = [0.0, 0.0, 0.0, 0.0]

  for i in 0..<values.len:
    result[0] += values[i]
  result[0] /= values.len.float

  for i in 0..<values.len:
    result[1] += (values[i] - result[0])^2
  result[1] /= values.len.float
  let sigma = sqrt(result[1])
  
  for i in 0..<values.len:
    let stdx = (values[i] - result[0]) / sigma
    result[2] += stdx^3
    result[3] += stdx^4
  result[2] /= values.len.float
  result[3] /= values.len.float


proc emitHeader(destination: File, output_format: OutputFormat) = 
  if output_format == OutputFormat.Complete:
    destination.write("model\tsample\titer\tseed\tnafbins\texternal_ref_af\tvalue\n")
  else:
    destination.write("model\titer\tseed\tnafbins\texternal_ref_af\tm1\tm2\tm3\tm4\n")


proc emitValues(model: Model, dosages: Dosages, iter: int, seed: int, n_afbins: int, values: seq[float], output_format: OutputFormat, af_source: SamplingReference, destination: File) = 
  let ext_ref_af = if af_source == SamplingReference.External: '1' else: '0'
  if output_format == OutputFormat.Complete:
    for i in 0..<values.len:
      destination.write(model.id & "\t" & dosages.samples[i] & "\t" & $iter & "\t" & $seed & "\t" & $n_afbins & "\t" & ext_ref_af & "\t" & $values[i] & "\n")
  else:
    let summary = calcSummaryStatistics(values)
    destination.write(model.id & "\t" & $iter & "\t" & $seed & "\t" & $n_afbins & "\t" & ext_ref_af & "\t" & $summary[0] & "\t" & $summary[1] & "\t" & $summary[2] & "\t" & $summary[3] & "\n")


proc calculationLoop(dosage_path: string, model_path: string, output_file: File, output_format: OutputFormat, af_source: SamplingReference, iters: int, n_afbins: int, seed: int) = 
  let dosages = loadDosages(dosage_path, n_afbins)
  let models = loadModels(model_path, n_afbins, dosages)

  stderr.write("Preparing random seed vector...\n")
  randomize(seed)
  var subseeds: seq[int] = @[]
  for i in 0..<iters:
    subseeds.add(random(int.high))

  stderr.write("Writing output...\n")
  emitHeader(output_file, output_format)

  var j = 0
  for model_id, model in models.pairs:
    j += 1
    stderr.eraseLine()
    stderr.write("Model " & $j & " / " & $models.len & ": " & model_id)
    let native_values = calcValues(model, dosages)
    emitValues(model, dosages, 0, seed, n_afbins, native_values, output_format, af_source, output_file)
    for i in 1..iters:
      let null_model = generateNullModel(model, dosages, af_source, subseeds[i-1])
      let null_values = calcValues(null_model, dosages)
      emitValues(model, dosages, i, seed, n_afbins, null_values, output_format, af_source, output_file)

  stderr.eraseLine()
  stderr.write("Done.\n")


proc printUsage(message: string = "") =
  if message != "":
    stderr.write(message & "\n\n")
  stderr.write("""
popstar: Calculate polygenic models and permuted nulls.

Usage: popstar [options] --dosages|d=DOSAGES --models|m=MODELS

Required parameters:
  --dosages|d=DOSAGES  Path to the input allele dosages file
  --models|m=MODELS    Path to the input model coefficient file

Options:
  --out|o=OUT     Path to the output file [default: stdout]
  --format|f=FMT  Output format (complete or summary) [default: complete]
  --iter|i=ITER   Number of resampled null iterations to calculate [default: 1000]
  --ref|r=REF     Source of resampling target allele frequency (external or internal) [default: external]
  --bins|b=BINS   Number of allele frequency bins for null allele matching [default: 100]
  --seed|s=SEED   PRNG seed [default: 314159265]

v""" & VERSION & "\n" & AUTHOR & "\n\n")


proc main() =
  var
    dosage_path: string = nil
    model_path: string = nil
    output_path: string = nil
    seed: int = 314159265
    iters: int = 1000
    bins: int = 100
    af_source: SamplingReference = SamplingReference.External
    output_format: OutputFormat = OutputFormat.Complete

  for kind, key, val in getopt():
    case kind
    of cmdArgument:
      printUsage("ERROR: does not accept positional arguments.")
      return
    of cmdLongOption, cmdShortOption:
      case key
      of "bins", "b": bins = val.parseInt
      of "dosages", "d": dosage_path = val
      of "format", "f":
        if val == "complete":
          output_format = OutputFormat.Complete
        elif val == "summary":
          output_format = OutputFormat.Summary
        else:
          printUsage("ERROR: parameter to --format|f must be either \"complete\" or \"summary\" (without quotes); was supplied: \"" & val & "\"")
          return
      of "iter", "i": iters = val.parseInt
      of "models", "m": model_path = val
      of "out", "o": output_path = val
      of "ref", "r":
        if val == "external":
          af_source = SamplingReference.External
        elif val == "internal":
          af_source = SamplingReference.Internal
        else:
          printUsage("ERROR: parameter to --ref|r must be either \"external\" or \"internal\" (without quotes); was supplied: \"" & val & "\"")
          return
      of "seed", "s": seed = val.parseInt
      else:
        printUsage("ERROR: unrecognised option " & key)
        return
    else: raise newException(Exception, "Unexpected kind " & $kind & " returned by getopt -- please report this error")

  if dosage_path == nil:
    printUsage("ERROR: Dosage file path is required.")
    return

  if model_path == nil:
    printUsage("ERROR: Model file path is required.")
    return

  let output_file = if output_path == nil: stdout else: system.open(output_path, fmWrite)

  calculationLoop(dosage_path, model_path, output_file, output_format, af_source, iters, bins, seed)


main()
