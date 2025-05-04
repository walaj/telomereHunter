# TelomereHunter (Modified Version)

This is a modified version (now with Python3 compatibility) of [TelomereHunter](https://doi.org/10.1186/s12859-019-2851-0), originally developed by Lars Feuerbach, Philip Ginsbach, and Lina Sieverling, and published in the paper:

**Feuerbach, L., Sieverling, L., Deeg, K.I. et al. TelomereHunter – in silico estimation of telomere content and composition from cancer genomes. BMC Bioinformatics 20, 272 (2019).**  
[https://doi.org/10.1186/s12859-019-2851-0](https://doi.org/10.1186/s12859-019-2851-0)

## Authors:
- Lars Feuerbach <l.feuerbach@dkfz-heidelberg.de>
- Philip Ginsbach
- Lina Sieverling <l.sieverling@dkfz-heidelberg.de>

## Modifications:
This repository contains modifications to the original TelomereHunter code. The changes made are listed below:
- Updated Pythxon2 -> Python3 syntax
- Used Python3 PySam language for reading BAM files

## Usage
Usage is provided with command:
```bash
PYTHONPATH=~/git python ~/git/telomereHunter/bin/telomereHunter -h
```

```bash
telomereHunter [-h] [-ibt TUMOR_BAM] [-ibc CONTROL_BAM] -o OUTPUT_DIR -p PID [-b BANDING_FILE] [-rt REPEAT_THRESHOLD] [-mqt MAPQ_THRESHOLD] [-d]
               [-r REPEATS [REPEATS ...]] [-con] [-gc1 LOWERGC] [-gc2 UPPERGC] [-nf] [-pl] [-pff {pdf,png,svg,all}] [-p1] [-p2] [-p3] [-p4] [-p5] [-p6]
               [-prc]
```	       

## License
This project is licensed under the GPL License. You can freely use, modify, and distribute this code under the terms of the GPL, as described in the LICENSE file.

## License Notice
This modified version of TelomereHunter retains the original GPL license as outlined in the [LICENSE](LICENSE) file. The code in this repository is provided as-is, with no warranty, express or implied.

Please see the LICENSE file for more details on the GPL license terms.

## Acknowledgements
This software was developed by Lars Feuerbach, Philip Ginsbach, and Lina Sieverling. The original work was supported by the research conducted at DKFZ Heidelberg.

For citation details, please refer to the original publication:
**Feuerbach, L., Sieverling, L., Deeg, K.I. et al. TelomereHunter – in silico estimation of telomere content and composition from cancer genomes. BMC Bioinformatics 20, 272 (2019).**  
[https://doi.org/10.1186/s12859-019-2851-0](https://doi.org/10.1186/s12859-019-2851-0)
