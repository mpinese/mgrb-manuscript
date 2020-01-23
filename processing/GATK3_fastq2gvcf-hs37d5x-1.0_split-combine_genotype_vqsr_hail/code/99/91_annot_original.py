#!./bin/pyhail.sh
import hail

hc = hail.HailContext(tmp_dir = 'tmp/hail')


def getAnnType(annotation, schema):
    ann_path = annotation.split(".")[1:]
    ann_type = schema
    for p in ann_path:
        ann_type = [x for x in ann_type.fields if x.name == p][0].typ
    return ann_type


def annotate_non_split_from_split(hc, non_split_vds_path, split_vds, annotations):
    ann_types = list(map(lambda x: str(getAnnType(x, split_vds.variant_schema)), annotations))

    variant_annotated_vds = (
        hc.read(non_split_vds_path, drop_samples=True)
        .annotate_variants_expr('va.variant = v')
        .split_multi()
    )

    ann_agg_codes = ["`%s` = index(va.map(x => {val: %s, aIndex: va.aIndex}).collect(), aIndex)" % (a, a) for a in annotations]
    agg = (
        split_vds
            .annotate_variants_vds(variant_annotated_vds, 'va.variant = vds.variant, va.aIndex = vds.aIndex')
            .filter_variants_expr('isDefined(va.variant)')
            .variants_table()
            .aggregate_by_key('variant = va.variant', ann_agg_codes + ['nAltAlleles = va.map(x => x.aIndex).max()'])
    )

    ann_codes = ['%s = let x = table.`%s` in' \
                 ' range(table.nAltAlleles).map(i => if(x.contains(i+1)) x[i+1].val else NA: %s)' % (a, a, b)
                 for (a, b) in zip(annotations, ann_types)]

    return (
        hc.read(non_split_vds_path)
            .annotate_variants_table(agg, expr=",".join(ann_codes))
    )



vds_annots = hc.read('test_annots.vds')


annotated_vds = annotate_non_split_from_split(hc, 'test_vars.vds', vds_annots, ['va.info.annot'])

annotated_vds.export_vcf('test_vars_annotorig.vcf.bgz')
annotated_vds.export_variants('test_vars_annotorig.tsv', 'variant = v, va.info.*')
