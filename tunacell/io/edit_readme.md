!Automatic readme editing file for data export.
!
!
# Lineage export README

Text file is structured with HEADER and CONTENT.

## HEADER

HEADER is marked with '!' as comment character. You can read informations such
as:
* experiment metadata
* kind of decomposition
    - `Independent lineages`: reconstructed trees are decomposed in independent
      lineages, such that each cell in the tree appears only once in the file.
    - `All lineages`: all paths from root to leaves for each reconstructed tree.
    - date of processing
* Filter set that have been applied to data, from reading from flat files to
  writing this output file. Note that filters applied in the image processing
  are not reported here, but may be found in the experiment metadata (if it has
  been given).


## A WORD ON DATA ORGANIZATION 

The organization of data in the associated file follows the hierarchical
pattern:

    EXPERIMENT > CONTAINER > TREE > LINEAGE > CELL

A given experiment is usually divided in smaller units, called CONTAINER.
It can be mapped onto a given experimental field of view, or a subset thereof
(depends on the experimental setup). It corresponds usually to the smaller
undivisible piece of image where segmentation and tracking have been performed.
Within such container, cell identifiers are unique; they may not be between two
or more containers.

## CONTENT

The division in container is reported in this file export by marking new
containers with a new line starting with keyword `CONTAINER`, preceded and
followed by empty lines. Example:

    CONTAINER: 21_0001

Within each container, few trees have been reconstructed. When a new tree is
reported, it starts with a new line with keyword `TREE`, preceded by one empty
line, and the identifier of root cell is given. Example:

    TREE (rootid: 3)

A given `TREE` is divided in independent lineages of cells. Each new lineage is
reported on a new line starting with keyword `LINEAGE`, preceded by one empty
line, and the sequence of cell identifiers is given.  Example:

    LINEAGE: 3,9,23

Lastly, each cell that belongs to the given independent lineage is reported by
a new line that starts with keyword `CELL` (no empty line separates cells from
same lineage).

    CELL:9, tau(mins):55.0, mu(db/hr):1.5496e+00, V_i(um^3):7.5167e-01,
    V_f(um^3):3.2847e+00

## CELL

Each `CELL` entry consists of a first line starting with keyword `CELL`,
followed by the cell identifier. Then, when possible, aggregate data is shown,
such as:
* `tau`: interdivision time (in minutes)
* `mu`: cell cycle estimate of volume growth rate in terms of doublings per
  hour. `mu` is obtained by performing a polynomial fit of order 1 of the log
  of volumes w.r.t times, divided by log(2.) (to get number of doublings) and
  multiplied by 60. (to get it in /hour).
* `V_i`: estimate of initial volume. Extrapolation of volume at mid-point
  between parent cell's last frame and cell's first frame, using a fit of log
  volume over 3 consecutive points (first 3 frames). In um^3.
* `V_f`: estimate of final volume. Extrapolation of volume at mid-point between
  cell's last frame and daugther cells' first frame, using a fit of log volume
  over 3 consecutive points (last 3 frames). In um^3.

Then each available observable is reported every given line by tab separated
values:
* time: minutes
* length: micrometres (um)
* width: um
* area: um^2 (projected area on horizontal plane)
* fluo: fluorescence in arbitrary units
* xCM: x position of Center of Mass, in um
* yCM: y position of Center of Mass, in um
* ALratio: area/length ratio, as an alternative report of lateral dimension
* volume: estimated from length and width assuming sphero-cylinder hypothesis
* concentration: fluo/volume
* density: fluo/area
* age: (also called phase) time spent in cell cycle from birth, divided by
  interdivision time. At birth age=0, at division age=1

Example:

    CELL: identifier
    time    time_frame0	time_frame1 time_frame2	...
    length  l_frame0	l_frame1    l_frame2	...
    ...

