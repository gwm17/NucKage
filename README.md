# NucKage

NucKage is a simulation of nuclear reaction kinematics and detector efficiency, with an interest towards simulating multiple reaction chains (i.e. a production reaction and then subsequent decays). Inspired by the work of one Nicolas Cage, this simulation aims to have a broad range (of applications) while still bringing (computing) performance.

>"Got my first chemistry set when I was seven, blew my eyebrows off, we never saw the cat again, been into it ever since."
>
>Nicolas Cage, *The Rock*

## Installation
Simply clone the repository using `git clone https://github.com/gwm17/NucKage.git`. NucKage uses the *[premake5](https://premake.github.io/)* build system. Download and install premake, and then generate the NucKage build files for your system (MacOS, Linux, Windows) using premake. NucKage has one external dependence, *[CERN's ROOT](https://root.cern.ch/)*. Since ROOT does not by default suport premake5 as a build system, the ROOT include directories and library directories must be manually entered into the NucKage premake5.lua file. Your ROOT directories can be found using the `root-config` script in MacOS and Linux. A typical Linux process might look like

- `git clone https://github.com/gwm17/NucKage.git`
- `cd NucKage`
- `premake5 gmake2`
- `make -j 4`

NucKage comes with a UI to generate configuration files, called Roles. The RoleGUI is written in python and uses Qt5 with the qtpy front-end wrapper. To use the RoleGUI one must have installed the qtpy library as well as one of the supported QT5 libraries (pyqt5 or PySide2). To launch the RoleGUI simply run `./bin/RoleGUI` from the top level directory of the repository.

## Usage
NucKage expects to be run from the top level directory of the repository as `./bin/NucKage <config>`. Currently NucKage accepts a single argument as the configuration (role) file to be used. Configurations are in plain-text, so with an example one would be able to write a role from scratch, however the RoleGUI is provided to make generating roles more straightforward as well as provide some simple checks to make sure a role will actually be valid for NucKage. In the future NucKage will also allow for an argument specifying the number of threads to allocate to the thread pool. NucKage saves a set of histograms and graphs to a ROOT outputfile specified in the configuration file.

## Principles

### Kinematics
NucKage is a flexible kinematics simulation. It is flexible in that it can handle a variety of reaction scenarios. For example, the reaction 10B(3He, a) populates states of 9B which can potentially decay to a system of p+8Be. That 8Be nucleus can then further decay to two alpha products. One would want to simulate this entire series of reactions, or in the language of NucKage this reaction chain. In general NucKage should be capable of handling any type of reaction chain, with some small caveats. First, NucKage assumes that the residual (i.e. final product) of the previous reaction in the chain serves as the target (or parent) of the following reaction. This is not flexible. Second, reactions involving gamma rays are not well defined by NucKage currently. NucKage will log errors if either of these rules are violated.

### Detector Efficiency
NucKage provides the framework to test detector geometric efficiencies. Detectors are described simply as geometries, and reaction products are tested for whether they fall within that geometry. In prinicple this framework could be expanded to include more physics based effects (energy loss, scattering, etc.). NucKage currently ships with the geometry for the SE-SPS aperature and SABRE array in use at FSU. Adding detector geomtry is relatively simple, and the included packages can be used as an example.

### Energy Loss
NucKage includes a built in energy loss integration framework. The code is based upon the work of Dale Visser, who wrote SPANC while working on the SE-SPS at Yale. The numerical methods are from Ziegler's SRIM (see the code for more details). The energy loss calculation is fairly expensive, and is the main bottleneck in the code. However, using an integration framework rather than a lookup table of values allows for NucKage to be capable of simulating energy loss through any solid material (within reason). These methods have been reported to be within approximately 5% accuracy, and this can be verified using the simple tests in the NucKage repository. There of course will be some cases for which these methods fail, so feel free to use the test to compare to other energy loss tools such as LISE++. NucKage currently does not have the ability to calculate energy loss for gases.

Energy loss is calculated for two kinds of particles: projectiles and ejectiles. All other particles (targets, residuals) are not used for energy loss. That is, in a chain like the 10B(3He, a) case mentioned above the 9B and 8Be residuals are not sent through energy loss while the 3He, alphas, and protons are. In essence, NucKage assumes that all reactions occur at the same location. The veracity of this assumption is up to the user to determine. 

### Performance
In general, nuclear physics experiments do not actually run a single reaction. A beam-like projectile is impinged upon a target and many possible reactions can take place. In order to properly understand the kinematics and detector performance, one would like to be able to run a simulation of all possible channels that are open in a uniform simulation environment. However, simulating so many reactions can be quite time consuming when running them one at a time (especially when striving to achieve an appropriate level of statistics).

In an effort to leverage modern hardware, NucKage utilizes a thread pool to run multiple simulations at the same time. In general, optimal performance will occur when there is a thread for every reaction, with gains as large as a factor of 2 in overal runtime; even in a worst case scenario (several reaction chains with many reactions using a single worker thread) performance gains are non-neglible as NucKage still utilizes the main thread to handle plotting of results while the worker thread runs the actual simulation. For insights on how the thread pool is implemented, see src/ThreadPool.h.

### Adding new detector geometries
In principle, any type of detector geometry can be programed into NucKage by following the examples given of the SPS aperature and the SABRE array. The difficulty arises in that currently each detector is distinct and needs to be added to the detector array independently. Additionally, the RoleGUI and configuration files are not terribly easy to modify.

## Notes
NucKage is still under rapid development and will continue to be so until further notice.