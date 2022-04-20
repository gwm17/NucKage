# NucKage

NucKage is a simulation of nuclear reaction kinematics, with an interest towards simulating multiple reaction chains (i.e. a production reaction and then subsequent decays). Inspired by the work of one Nicolas Cage, this simulation aims to have a broad range (of applications) while still bringing (computing) performance.

>"Got my first chemistry set when I was seven, blew my eyebrows off, we never saw the cat again, been into it ever since."
>Nicolas Cage, *The Rock*

## Installation
Simply clone the repository using `git clone https://github.com/gwm17/NucKage.git`. NucKage uses the *[premake5](https://premake.github.io/)* build system. Download and install premake, and then generate the NucKage build files for your system (MacOS, Linux, Windows) using premake. NucKage has one external dependence, *[CERN's ROOT](https://root.cern.ch/)*. Since ROOT does not by default suport premake5 as a build system, the ROOT include directories and library directories must be manually entered into the NucKage premake5.lua file. Your ROOT directories can be found using the `root-config` script in MacOS and Linux. A typical Linux process might look like

- `git clone https://github.com/gwm17/NucKage.git`
- `cd NucKage`
- `premake5 gmake2`
- `make -j 4`

## Usage
NucKage expects to be run from the top level directory of the repository (as `./bin/NucKage <config>`). Currently NucKage accepts a single argument as the configuration file to be used. Configurations are currently written by hand, however a UI is planned to be developed to generate config files. 