# Project: Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary Ulva sp. supplementation: A meta-analysis 

## Step2: Example script for extracting data from figures in original articles. 

## Relevant figures were saved as image files and are stiored in "Data/figures_for_extraction" folder.
## They are loaded one by one via MetaDigitise package for manual extraction of relevant values.
## Output is stored as R object for incorporation in teh main data set.

## Install packages
Install.packages(c("MetaDigitise", "here"))

## Load required libraries
library(here)
library(MetaDigitise)


## Load images from first paper.

image_dir_Falade <- here("UlvaSupplements_AbaloneMA", "Data", "figures_for_extraction", "Falade")

## Start  digitising process

### Note: this is done ion the Console and Plots window of RStudio, not as executable code.
### The following lines are the inputs provided during the digitising process, for the records.


digitised_data_Falade <- metaDigitise(image_dir_Falade, summary = TRUE)
print(digitised_data_Falade)

#Image processing
#S006.1_CF
1
d
c
m
final condition factor
0.98
1.14
n
y
control
5
a
IMTA
5
a
non-IMTA
5
f
CI95
n

#S006.1_FinalWeight
1
d
c
m
Final weight (g)
100
122
n
y
control
5
a
IMTA
5
a
non-IMTA
5
f
CI95
n

#S006.1_SGR
1
d
c
m
SGR
0.17
0.26
n
y
control
5
a
IMTA
5
a
non-IMTA
5
f
CI95
n

#S006.2_FCR
y
c
m
FCR
1.0
4.5
n
y
IMTA-U0
3
a
IMTA-U0.75
3
a
IMTA-U1.5
3
a
IMTA-U3
3
a
IMTA-U6
3
a
IMTA-U12
3
f
CI95
y

#S006.2_FEED INTAKE
c
m
FEED INTAKE
5.0
11.0
n
y
IMTA-U0
3
a
IMTA-U0.75
3
a
IMTA-U1.5
3
a
IMTA-U3
3
a
IMTA-U6
3
a
IMTA-U12
3
f
CI95
y

#S006.2_PER
c
m
PER
1.0
3.0
n
y
IMTA-U0
3
a
IMTA-U0.75
3
a
IMTA-U1.5
3
a
IMTA-U3
3
a
IMTA-U6
3
a
IMTA-U12
3
f
CI95
y

#Image proccessing
#S006.2_CF
1
s
s
c
CF
Ulva inclusion
0
1.6
0
15
n
n
IMTA-U0
IMTA-U0.75
IMTA-U1.5
IMTA-U3
IMTA-U6
IMTA-U12
c
a
f
n
y

#Image proccessing
#S006.2_FINAL SHELL LENGTH
c
FSL
Ulva inclusion
68
76
0
15
n
n
IMTA-U0
IMTA-U0.75
IMTA-U1.5
IMTA-U3
IMTA-U6
IMTA-U12
c
a
f
n
y

#Image proccessing
#S006.2_SGR
c
SGR
Ulva inclusion
0
0.40
0
15
n
n
IMTA-U0
IMTA-U0.75
IMTA-U1.5
IMTA-U3
IMTA-U6
IMTA-U12
c
a
f
n
y

# Please download directory and run the following
here()
image_dir_Mwu <- here("UlvaSupplements_AbaloneMA", "Data", "figures_for_extraction", "Mwangudza")

# Start the digitising interface
digitised_data_Mwu <- metaDigitise(image_dir_Mwu, summary = TRUE)
print(digitised_data_Mwu)

#Final length
1
d
c
m
length(mm)
54
70
n
y
Abfeed (Reference) 30 days
4
a
Abfeed (Reference) 120 days
4
a
Abfeed (Reference) 150 days
4
a
Biosecure Ulva 30 days
4
a
Biosecure Ulva 120 days
4
a
Biosecure Ulva 150 days
4
a
Not biosecure Ulva 30 days
4
a
Not biosecure Ulva 120 days
4
a
Not biosecure Ulva 150 days
4
f
se
n

#Final weight
1
d
c
m
weight(g)
0
65
n
y
Abfeed (Reference) 30 days
4
a
Abfeed (Reference) 120 days
4
a
Abfeed (Reference) 150 days
4
a
Biosecure Ulva 30 days
4
a
Biosecure Ulva 120 days
4
a
Biosecure Ulva 150 days
4
a
Not biosecure Ulva 30 days
4
a
Not biosecure Ulva 120 days
4
a
Not biosecure Ulva 150 days
4
f
se
n
y
