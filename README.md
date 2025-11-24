# LoPhy

**LoPhy** is an algorithm designed for **LO**ngitudinal **Phy**logenetic reconstruction from single-cell amplicon sequencing. The algorithmic design and experimental validation can be found here: [https://www.biorxiv.org/content/10.1101/2025.09.16.676596v1.abstract](https://www.biorxiv.org/content/10.1101/2025.09.16.676596v1.abstract "https://www.biorxiv.org/content/10.1101/2025.09.16.676596v1.abstract").

### Requirements

* MacOS or Linux
* \>= c++20 

### Installing
`git clone https://github.com/ethanumn/LoPhy.git`
`cd LoPhy`
`mkdir bin`
`make`


### Run on example data

Using this command will produce the 
```
cd LoPhy
./bin/LoPhy -c example/inputs/character_matrix.csv -v example/inputs/variant_reads.csv -t example/inputs/total_reads.csv -m example/inputs/meta.csv -r example/inputs/region_reads.csv -o example/outputs -p out -seed 0 -homp 15 -hetp 4 -fp 0.02 -fn 0.05 -s example/inputs/cell_samples.txt
```
