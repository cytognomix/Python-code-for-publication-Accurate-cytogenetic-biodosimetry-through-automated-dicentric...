# Python-code-for-publication-Accurate-cytogenetic-biodosimetry-through-automated-dicentric...
Image Selection for Metaphase Images in Automated Dicentric Chromosome Identification
Author: Yanxin Li, Jin Liu

List of files and folders:
DemoSample1: metaphase images for DemoSample1
DemoSample2: metaphase images for DemoSample2
DemoSample1.adcisample: ADCI sample file generated by ADCI software for DemoSample1.
DemoSample2.adcisample: ADCI sample file generated by ADCI software for DemoSample2.
ADCI sample files contain area and filters information of images which is required by Image Selection code.
DemoScript.py: Image Selection code

Prerequisites for the Image Selection code:
Python 2. Code is tested in Python 2.7.6
NumPy

How to run:
Run DemoScript.py. Make sure DemoSample1.adcisample and DemoSample2.adcisample are in the same folder as DemoScript.py.

Results:
For each sample, the code will display image names in descending order of quality, measured by Group Bin Distance or Combined Z Score (using weight 434521).
To evaluate the sorting, view metaphase images in DemoSample1 or DemoSample2.