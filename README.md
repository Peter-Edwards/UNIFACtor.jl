# UNIFACtor.jl

UNIFACtor.jl is a package for Julia that allows the easy calculation of activity coefficients of an arbitrary number of different components in a liquid mixture.

The current iteration of the package includes the original UNIFAC method and the Dortmund Modified method to find the activity coefficients. Any and all of the UNIFAC coefficients used in the model were taken from the official list published by the UNIFAC Consortium:

http://unifac.ddbst.de/published-parameters-unifac.html
http://unifac.ddbst.de/PublishedParametersUNIFACDO.html

the first 50 groups UNIFAC groups are included in the package’s UNIFAC model along with 108 different sub groups.
All 63 groups are included in the modified UNIFAC model, along with the requisite 125 sub groups

## Usage

Currently, there are 2 functions avaliable for use:

```julia
UNIFAC.Activity(T_k,M_lst,x_arr)
UNIFACmod.Activity(T_k,M_lst,x_arr)
```

T_k is a temperature value in kelvin. The x_arr is a 1xn matrix of molar faction values (where n is the number of components). M_lst is a vector of 1xm matrices describing the composition of each of the component molecules, where m is dependant on the functional groups present. As stated previously, 108 different functional groups are available for use, however the matrices present in in M_lst only need to be long enough to fully describe the molecule. An example of how to use the UNIFAC function is given below:

```julia
x1=0.5
x2=0.25
x3=1-x1-x2
x_arr=[x1 x2 x3]

T_k=298

M_1=[0 0 0 0 0 0 0 0 0 6]; #benzene
M_2=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1]; #ethanol
M_3=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone

M_lst=[M_1,M_2,M_3]

UNIFAC.Activity(T_k,M_lst,x_arr)
```
In each case it is assumed that all of the remaining functions groups out of the 108 are treated as zero. The modified UNIFAC function also has this feature. Please note that some molecule composition matrices may need to be changed if switching from UNIFAC to modified UNIFAC and vice versa because each method has their own set of groups and subgroups, for example:

```julia
M_3=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone in UNIFAC
M_3=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone in UNIFACmod
```
## Data Tables
The order in which the functional groups are listed for the input matrices is different to the order listed of the official list of parameters. The order in which they are arranged in both the UNIFAC and modified UNIFAC is based solely on their interaction groups. The List of Usable subgroups and groups is listed below. For creating the input molecule composition matrices, use the UNIFACtor Index column.


### Basic UNIFAC function:

| UNIFACtor Index | Original Index | Sub Group |    Group     |
|-----------------|----------------|-----------|--------------|
|               1 |              1 | CH3       | [1]CH2       |
|               2 |              2 | CH2       | [1]CH2       |
|               3 |              3 | CH        | [1]CH2       |
|               4 |              4 | C         | [1]CH2       |
|               5 |              5 | CH2=CH    | [2]C=C       |
|               6 |              6 | CH=CH     | [2]C=C       |
|               7 |              7 | CH2=C     | [2]C=C       |
|               8 |              8 | CH=C      | [2]C=C       |
|               9 |             70 | C=C       | [2]C=C       |
|              10 |              9 | ACH       | [3]ACH       |
|              11 |             10 | AC        | [3]ACH       |
|              12 |             11 | ACCH3     | [4]ACCH2     |
|              13 |             12 | ACCH2     | [4]ACCH2     |
|              14 |             13 | ACCH      | [4]ACCH2     |
|              15 |             14 | OH        | [5]OH        |
|              16 |             15 | CH3OH     | [6]CH3OH     |
|              17 |             16 | H2O       | [7]H2O       |
|              18 |             17 | ACOH      | [8]ACOH      |
|              19 |             18 | CH3CO     | [9]CH2CO     |
|              20 |             19 | CH2CO     | [9]CH2CO     |
|              21 |             20 | CHO       | [10]CHO      |
|              22 |             21 | CH3COO    | [11]CCOO     |
|              23 |             22 | CH2COO    | [11]CCOO     |
|              24 |             23 | HCOO      | [12]HCOO     |
|              25 |             24 | CH3O      | [13]CH2O     |
|              26 |             25 | CH2O      | [13]CH2O     |
|              27 |             26 | CHO       | [13]CH2O     |
|              28 |             27 | THF       | [13]CH2O     |
|              29 |             28 | CH3NH2    | [14]CNH2     |
|              30 |             29 | CH2NH2    | [14]CNH2     |
|              31 |             30 | CHNH2     | [14]CNH2     |
|              32 |             31 | CH3NH     | [15]CNH      |
|              33 |             32 | CH2NH     | [15]CNH      |
|              34 |             33 | CHNH      | [15]CNH      |
|              35 |             34 | CH3N      | [16]\(C\)3N    |
|              36 |             35 | CH2N      | [16]\(C\)3N    |
|              37 |             36 | ACNH2     | [17]ACNH2    |
|              38 |             37 | C5H5N     | [18]PYRIDINE |
|              39 |             38 | C5H4N     | [18]PYRIDINE |
|              40 |             39 | C5H3N     | [18]PYRIDINE |
|              41 |             40 | CH3CN     | [19]CCN      |
|              42 |             41 | CH2CN     | [19]CCN      |
|              43 |             42 | COOH      | [20]COOH     |
|              44 |             43 | HCOOH     | [20]COOH     |
|              45 |             44 | CH2CL     | [21]CCL      |
|              46 |             45 | CHCL      | [21]CCL      |
|              47 |             46 | CCL       | [21]CCL      |
|              48 |             47 | CH2CL2    | [22]CCL2     |
|              49 |             48 | CHCL2     | [22]CCL2     |
|              50 |             49 | CCL2      | [22]CCL2     |
|              51 |             50 | CHCL3     | [23]CCL3     |
|              52 |             51 | CCL3      | [23]CCL3     |
|              53 |             52 | CCL4      | [24]CCL4     |
|              54 |             53 | ACCL      | [25]ACCL     |
|              55 |             54 | CH3NO2    | [26]CNO2     |
|              56 |             55 | CH2NO2    | [26]CNO2     |
|              57 |             56 | CHNO2     | [26]CNO2     |
|              58 |             57 | ACNO2     | [27]ACNO2    |
|              59 |             58 | CS2       | [28]CS2      |
|              60 |             59 | CH3SH     | [29]CH3SH    |
|              61 |             60 | CH2SH     | [29]CH3SH    |
|              62 |             61 | FURFURAL  | [30]FURFURAL |
|              63 |             62 | DOH       | [31]DOH      |
|              64 |             63 | I         | [32]I        |
|              65 |             64 | BR        | [33]BR       |
|              66 |             65 | CH=-C     | [34]C=-C     |
|              67 |             66 | C=-C      | [34]C=-C     |
|              68 |             67 | DMSO      | [35]DMSO     |
|              69 |             68 | ACRY      | [36]ACRY     |
|              70 |             69 | CL-(C=C)  | [37]CLCC     |
|              71 |             71 | ACF       | [38]ACF      |
|              72 |             72 | DMF       | [39]DMF      |
|              73 |             73 | HCON(..   | [39]DMF      |
|              74 |             74 | CF3       | [40]CF2      |
|              75 |             75 | CF2       | [40]CF2      |
|              76 |             76 | CF        | [40]CF2      |
|              77 |             77 | COO       | [41]COO      |
|              78 |             78 | SIH3      | [42]SIH2     |
|              79 |             79 | SIH2      | [42]SIH2     |
|              80 |             80 | SIH       | [42]SIH2     |
|              81 |             81 | SI        | [42]SIH2     |
|              82 |             82 | SIH2O     | [43]SIO      |
|              83 |             83 | SIHO      | [43]SIO      |
|              84 |             84 | SIO       | [43]SIO      |
|              85 |             85 | NMP       | [44]NMP      |
|              86 |             86 | CCL3F     | [45]CCLF     |
|              87 |             87 | CCL2F     | [45]CCLF     |
|              88 |             88 | HCCL2F    | [45]CCLF     |
|              89 |             89 | HCCLF     | [45]CCLF     |
|              90 |             90 | CCLF2     | [45]CCLF     |
|              91 |             91 | HCCLF2    | [45]CCLF     |
|              92 |             92 | CCLF3     | [45]CCLF     |
|              93 |             93 | CCL2F2    | [45]CCLF     |
|              94 |             94 | AMH2      | [46]CON(AM)  |
|              95 |             95 | AMHCH3    | [46]CON(AM)  |
|              96 |             96 | AMHCH2    | [46]CON(AM)  |
|              97 |             97 | AM(CH3)2  | [46]CON(AM)  |
|              98 |             98 | AMCH3CH2  | [46]CON(AM)  |
|              99 |             99 | AM(CH2)2  | [46]CON(AM)  |
|             100 |            100 | C2H5O2    | [47]OCCOH    |
|             101 |            101 | C2H4O2    | [47]OCCOH    |
|             102 |            102 | CH3S      | [48]CH2S     |
|             103 |            103 | CH2S      | [48]CH2S     |
|             104 |            104 | CHS       | [48]CH2S     |
|             105 |            105 | MORPH     | [49]MORPH    |
|             106 |            106 | C4H4S     | [50]THIOPHEN |
|             107 |            107 | C4H3S     | [50]THIOPHEN |
|             108 |            108 | C4H2S     | [50]THIOPHEN |

### Dortmund Modified UNIFAC function:

| UNIFACtor Index | Original Index | Subgroup |    Group     |
|-----------------|----------------|----------|--------------|
|               1 |              1 | CH3      | [1]CH2       |
|               2 |              2 | CH2      | [1]CH2       |
|               3 |              3 | CH       | [1]CH2       |
|               4 |              4 | C        | [1]CH2       |
|               5 |              5 | CH2=CH   | [2]C=C       |
|               6 |              6 | CH=CH    | [2]C=C       |
|               7 |              7 | CH2=C    | [2]C=C       |
|               8 |              8 | CH=C     | [2]C=C       |
|               9 |             70 | C=C      | [2]C=C       |
|              10 |              9 | ACH      | [3]ACH       |
|              11 |             10 | AC       | [3]ACH       |
|              12 |             11 | ACCH3    | [4]ACCH2     |
|              13 |             12 | ACCH2    | [4]ACCH2     |
|              14 |             13 | ACCH     | [4]ACCH2     |
|              15 |             14 | OH (P)   | [5]OH        |
|              16 |             81 | OH (S)   | [5]OH        |
|              17 |             82 | OH (T)   | [5]OH        |
|              18 |             15 | CH3OH    | [6]CH3OH     |
|              19 |             16 | H2O      | [7]H2O       |
|              20 |             17 | ACOH     | [8]ACOH      |
|              21 |             18 | CH3CO    | [9]CH2CO     |
|              22 |             19 | CH2CO    | [9]CH2CO     |
|              23 |             20 | CHO      | [10]CHO      |
|              24 |             21 | CH3COO   | [11]CCOO     |
|              25 |             22 | CH2COO   | [11]CCOO     |
|              26 |             23 | HCOO     | [12]HCOO     |
|              27 |             24 | CH3O     | [13]CH2O     |
|              28 |             25 | CH2O     | [13]CH2O     |
|              29 |             26 | CHO      | [13]CH2O     |
|              30 |             28 | CH3NH2   | [14]CH2NH2   |
|              31 |             29 | CH2NH2   | [14]CH2NH2   |
|              32 |             30 | CHNH2    | [14]CH2NH2   |
|              33 |             85 | CNH2     | [14]CH2NH2   |
|              34 |             31 | CH3NH    | [15]CH2NH    |
|              35 |             32 | CH2NH    | [15]CH2NH    |
|              36 |             33 | CHNH     | [15]CH2NH    |
|              37 |             34 | CH3N     | [16]\(C\)3N    |
|              38 |             35 | CH2N     | [16]\(C\)3N    |
|              39 |             36 | ACNH2    | [17]ACNH2    |
|              40 |             37 | AC2H2N   | [18]PYRIDINE |
|              41 |             38 | AC2HN    | [18]PYRIDINE |
|              42 |             39 | AC2N     | [18]PYRIDINE |
|              43 |             40 | CH3CN    | [19]CH2CN    |
|              44 |             41 | CH2CN    | [19]CH2CN    |
|              45 |             42 | COOH     | [20]COOH     |
|              46 |             44 | CH2CL    | [21]CCL      |
|              47 |             45 | CHCL     | [21]CCL      |
|              48 |             46 | CCL      | [21]CCL      |
|              49 |             47 | CH2CL2   | [22]CCL2     |
|              50 |             48 | CHCL2    | [22]CCL2     |
|              51 |             49 | CCL2     | [22]CCL2     |
|              52 |             51 | CCL3     | [23]CCL3     |
|              53 |             52 | CCL4     | [24]CCL4     |
|              54 |             53 | ACCL     | [25]ACCL     |
|              55 |             54 | CH3NO2   | [26]CNO2     |
|              56 |             55 | CH2NO2   | [26]CNO2     |
|              57 |             56 | CHNO2    | [26]CNO2     |
|              58 |             57 | ACNO2    | [27]ACNO2    |
|              59 |             58 | CS2      | [28]CS2      |
|              60 |             59 | CH3SH    | [29]CH3SH    |
|              61 |             60 | CH2SH    | [29]CH3SH    |
|              62 |             61 | FURFURAL | [30]FURFURAL |
|              63 |             62 | DOH      | [31]DOH      |
|              64 |             63 | I        | [32]I        |
|              65 |             64 | BR       | [33]BR       |
|              66 |             65 | CH=-C    | [34]C=-C     |
|              67 |             66 | C=-C     | [34]C=-C     |
|              68 |             67 | DMSO     | [35]DMSO     |
|              69 |             68 | ACRY     | [36]ACRY     |
|              70 |             69 | CL-(C=C) | [37]CLCC     |
|              71 |             71 | ACF      | [38]ACF      |
|              72 |             72 | DMF      | [39]DMF      |
|              73 |             73 | HCON(..  | [39]DMF      |
|              74 |             74 | CF3      | [40]CF2      |
|              75 |             75 | CF2      | [40]CF2      |
|              76 |             76 | CF       | [40]CF2      |
|              77 |             77 | COO      | [41]COO      |
|              78 |             78 | CY-CH2   | [42]CY-CH2   |
|              79 |             79 | CY-CH    | [42]CY-CH2   |
|              80 |             80 | CY-C     | [42]CY-CH2   |
|              81 |             27 | THF      | [43]CY-CH2O  |
|              82 |             83 | CY-CH2O  | [43]CY-CH2O  |
|              83 |             84 | TRIOXAN  | [43]CY-CH2O  |
|              84 |             43 | HCOOH    | [44]HCOOH    |
|              85 |             50 | CHCL3    | [45]CHCL3    |
|              86 |             86 | NMP      | [46]CY-CONC  |
|              87 |             87 | NEP      | [46]CY-CONC  |
|              88 |             88 | NIPP     | [46]CY-CONC  |
|              89 |             89 | NTBP     | [46]CY-CONC  |
|              90 |             91 | CONH2    | [47]CONR     |
|              91 |             92 | CONHCH3  | [47]CONR     |
|              92 |            100 | CONHCH2  | [47]CONR     |
|              93 |            101 | AM(CH3)2 | [48]CONR2    |
|              94 |            102 | AMCH3CH2 | [48]CONR2    |
|              95 |            103 | AM(CH2)2 | [48]CONR2    |
|              96 |             93 | HCONHCH3 | [49]HCONR    |
|              97 |             94 | HCONHCH2 | [49]HCONR    |
|              98 |            104 | AC2H2S   | [52]ACS      |
|              99 |            105 | AC2HS    | [52]ACS      |
|             100 |            106 | AC2S     | [52]ACS      |
|             101 |            107 | H2COCH   | [53]EPOXIDES |
|             102 |            108 | COCH     | [53]EPOXIDES |
|             103 |            109 | HCOCH    | [53]EPOXIDES |
|             104 |            119 | H2COCH2  | [53]EPOXIDES |
|             105 |            153 | H2COC    | [53]EPOXIDES |
|             106 |            112 | (CH3)2CB | [55]CARBONAT |
|             107 |            113 | (CH2)2CB | [55]CARBONAT |
|             108 |            114 | CH2CH3CB | [55]CARBONAT |
|             109 |            110 | (CH2)2SU | [56]SULFONE  |
|             110 |            111 | CH2SUCH  | [56]SULFONE  |
|             111 |            122 | CH3S     | [61]CH2S     |
|             112 |            123 | CH2S     | [61]CH2S     |
|             113 |            124 | CHS      | [61]CH2S     |
|             114 |            178 | C3H2N2+  | [84]IMIDAZOL |
|             115 |            184 | C3H3N2+  | [84]IMIDAZOL |
|             116 |            179 | BTI-     | [85]BTI      |
|             117 |            189 | C4H8N+   | [87]PYRROL   |
|             118 |            195 | BF4-     | [89]BF4      |
|             119 |            196 | C5H5N+   | [90]PYRIDIN  |
|             120 |            220 | C5H4N+   | [90]PYRIDIN  |
|             121 |            197 | OTF-     | [91]OTF      |
|             122 |            201 | -S-S-    | [93]-S-S-    |
|             123 |            209 | SO4      | [98]SO4      |
|             124 |            210 | HSO4     | [98]SO4      |
|             125 |            211 | PF6      | [99]PF6      |
