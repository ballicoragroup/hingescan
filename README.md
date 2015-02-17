# hingescan

Hingescan is program for detecting hinges in protein structures.
Copyright 2015 Miguel A. Ballicora

## Compilation and Installation in GNU/Linux

type

`make`

to generate the program and then

`make install`

or `sudo make install` in Ubuntu. The program is very portable, so it should be very easy to compile in other systems.

## Usage

As an example, a typical use would be

`hingescan -a file1.pdb -b file2.pdb -w 3 -W 25 -o output.txt -g gnuplotfile.plt`

In this case, `file1.pdb` and `file2.pdb` will be compared with windows from 3 to 25 and the output, in matrix form will be in `output.txt`
To plot this information in a "heat map", anothe file will be saved (gnuplotfile.plt) that could be the input to be used with gnuplot.
Then, after the previous command was executed, the following (gnuplot needs to be already installed):

`gnuplot gnuplotfile.plt`

will yield the file `output.txt.svg`, which can be converted in linux to other graph formats. For instance with

`convert output.txt.svg output.jpg`

Another alternative is to produce hingescan scores with a give window. For instance

`hingescan -a file1.pdb -b file2.pdb -w 9 -o out.csv`

will output a comma separated value file (`out.csv`) that could be open with any spreadsheet program. The window size is 9.

There is a further help if hingescan is executed with the switch -h in the command line. That is:
```
Program to detect hinges by comparing two conformers

quick example 1: hingescan -a file1.pdb -b file2.pdb -w 21 -o out.csv
  - Processes file1.pdb and file2.pdb with a window size of 21 residues
  - Outputs a list of hinge scores for each residue in out.csv, which is
    a comma separated value file

quick example 2: hingescan -a file1.pdb -b file2.pdb -w 3 -W 21 -o out.txt -g gp.plt
  - Processes file1.pdb and file2.pdb from window 3 to 21 residues
  - Outputs a matrix of hinge scores for each residue in out.txt (text file)
  - Outputs a gnuplot 4.2 file (gp.plt) to plot out.txt in two dimensions.
    (score v. window). This file would be the input for gnuplot for a 2D plot

usage: hingescan [-OPTION]
 -h          print this help
 -H          print just the switches
 -v          print version number and exit
 -L          display the license information
 -q          quiet mode (no screen progress updates)
 -w <n>      window size (always odd). Min size if -W is provided
 -W <n>      multi scan, from -w to -W window sizes (always odd numbers)
 -Z <n>      top hinge score for plotting, default is the maximum observed
 -c          corrected hinge scores
 -a <file>   first input file in pdb format
 -b <file>   second input file in pdb format
 -o <file>   output file (text format), goes to the screen if not present
 -g <file>   output gnuplot 4.2 format, if -W switch is provided

Copyright (c) 2014 Miguel A. Ballicora
There is NO WARRANTY of any kind
```


