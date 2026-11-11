# Third-party licenses

Libraries compiled into `sage.exe` and `sage-gui.exe` whose licenses ask for their text and notices to travel with the programs:

* `sage-plus-raw/`: [sage-plus-raw](https://github.com/pgarrett-scripps/sage-plus-raw) 0.1.0, the Thermo RAW reader (Apache-2.0; Sage Plus's fork of [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw), Sigilweaver Holdings LLC). `LICENSE` and `NOTICE` are copied unchanged from the crate.

* `sage-plus/`: [Sage Plus](https://github.com/pgarrett-scripps/sage-plus) (MIT, the same license and copyright as Sage), from which XNova adapts its quality control, timsTOF calibration, nonlinear retention-time alignment, picked FDR, long-form LFQ table, named modifications, mass-offset modifications, site localization, site libraries and streamed prefilter. `LICENSE` is copied unchanged.

Methods and formats XNova reimplements without including their code:

* [EasyPQP](https://github.com/grosenberger/easypqp) (BSD-3-Clause, George Rosenberger and contributors): the method and columns of the spectral libraries (`spectral_library`).
* [xiDIA-library](https://github.com/Rappsilber-Laboratory/xiDIA-library) (Apache-2.0, Rappsilber Laboratory): the form of the crosslink spectral libraries (`xlibrary.*`) for Spectronaut.

Unimod's data and license are in `crates/sage-gui/assets/unimod/`; xiVIEW's are in the `xiview` folder of each release.
