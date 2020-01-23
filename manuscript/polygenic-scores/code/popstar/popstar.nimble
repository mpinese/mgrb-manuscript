# Package

version       = "0.0.2"
author        = "Mark Pinese <m.pinese@garvan.org.au>"
description   = "Poly/Oligogenic Phenotype Score Testing And Resampling"
license       = "MIT"


# Dependencies

requires "nim >= 0.17.0"


# Building

bin = @["popstar"]
skipExt = @["nim"]
skipDirs = @["tests", "utils"]


# Tests

task test, "Run the test suite":
  exec "nim c -r tests/tests"
