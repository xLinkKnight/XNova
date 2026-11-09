# Unimod, bundled with XNova

`unimod.xml` is [Unimod](https://www.unimod.org)'s XML export, unchanged, as
published at <https://www.unimod.org/xml/unimod.xml>:

| | |
|---|---|
| Last modified on unimod.org | 17 February 2026 |
| Downloaded | 9 November 2026 |
| Size | 2,506,678 bytes, 1,560 modifications (accessions up to #2147) |
| SHA-256 | `aa50802955726f15312773575c36324fd13ef82a1b0f9aa0e47015a4005c4e8d` |

`sage-gui` compiles this file in. Looking up a Unimod accession on the
Modifications tab therefore needs no internet connection and sends no traffic
to unimod.org. `sage-gui` only contacts unimod.org when you click **Check
unimod.org** for an accession that is not in this copy, and then at most once
a day. That download is kept in `%APPDATA%\sage-gui\unimod.xml` and is used
before this copy when it is newer.

## Copyright and license

Unimod is a public domain database, distributed under a copyleft licence, the
Design Science License. `LICENSE-unimod.txt` is that license, as published at
<https://www.unimod.org/dsl.txt>. The file's own notice reads:

> Copyright (C) 2002-2006 Unimod; this information may be copied, distributed
> and/or modified under certain conditions, but it comes WITHOUT ANY WARRANTY;
> see the accompanying Design Science License for more details

* XNova redistributes `unimod.xml` verbatim, with its copyright notice and this
  license (section 3 of the license). The `-text` entry in `.gitattributes`
  keeps git from changing its line endings.
* `sage-gui.exe` carries the same file compiled in (the license's "Object
  Form"). Its Source Data, this folder, is in XNova's repository and in each
  release download (section 3 (a)).
* XNova does not modify the file. The Design Science License covers Unimod
  only: aggregating it with XNova does not bring XNova's own code (MIT) under
  it (section 3).

## Citation

Creasy, D. M.; Cottrell, J. S. Unimod: Protein modifications for mass
spectrometry. *Proteomics* **2004**, 4 (6), 1534–1536.
<https://doi.org/10.1002/pmic.200300744>

## Updating

1. Download <https://www.unimod.org/xml/unimod.xml> once, and replace `unimod.xml` with it.
2. Update `BUNDLED_DATE` and `BUNDLED_MODIFIED` in `crates/sage-gui/src/unimod.rs`
   (unimod.org's `Last-Modified` for the file).
3. Update the table above.
