# RapidSim [![DOI](https://zenodo.org/badge/50026095.svg)](https://zenodo.org/badge/latestdoi/50026095)

RapidSim is a fast Monte Carlo generator for particle physics.
It can be useful for generating background, or even toy data sets
for acceptance studies.

* It uses [`TGenPhaseSpace`][tgenphasespace] to generate b events.
* [FONLL][fonll] is used to boost the b.
* The daughter particle momenta are smeared correctly.
* Cuts can be made on daughter particle properties.
* There are tools to swap mass hypothesis.


## Setup

Now moved to using CMake to compile a binary. To build do the following:

```shell
$ # Setup your ROOT environment. On lxplus, you can do this:
$ source /afs/cern.ch/sw/lcg/releases/LCG_85b/ROOT/6.06.06/x86_64-slc6-gcc49-opt/bin/thisroot.sh
$ source /afs/cern.ch/sw/lcg/releases/LCG_85b/GSL/2.1/x86_64-slc6-gcc49-opt/GSL-env.sh # required for ROOT (not RapidSim)
$ # Get the compiler (also works with clang). Again, on lxplus do this:
$ export CC=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_79/gcc/4.9.3/x86_64-slc6/bin/gcc
$ export CXX=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_79/gcc/4.9.3/x86_64-slc6/bin/g++
$ mkdir build
$ cd build
$ cmake ..
$ make -j 4
```

The usage is

```shell
$ src/RapidSim.exe <decay mode> <events to generate> <save tree?>
```

To run an example try

```shell
$ src/RapidSim.exe ../src/Bs2Jpsiphi 10000 1
```

To run the full system validation:

```shell
$ cd RapidSim/validation
$ source runValidation.sh
```

## Decays

To generate a new decay mode you must write a `.decay` file using the following syntax:

* particles are named according to `config/particles.dat` and separated by 
  spaces
* a decay is denoted by `->`
* subdecays are demarcated by braces `{}`, e.g.

```
Bs0 -> { Jpsi -> mu+ mu- } { phi -> K+ K- }
```

When RapidSim is first run for a new decay mode it will generate a simple `.config` file.
These files allow properties of the particles and global settings to be configured for each decay.
Run RapidSim with 0 events to generate if you just want to produce the `.config` file.
The default values listed below are those obtained when using the automatically generated `.config` file. 
Note that when lines are removed from the `.config` file this behaviour may change.

## Configuration

Global settings should be defined at the start of the file using the syntax:

```
<setting> : <value>
``````

Particle settings should be defined after the corresponding `@#` tag using the same syntax.

### Global settings

* `seed`:
  * Sets the seed for the random number generator
  * Default: 4357

* `acceptance`:
  * Sets the type of geometric acceptance to apply to the decay
  * Options are `Any`, `ParentIn`, `AllIn` and `AllDownstream`
  * Default: `Any`

* `geometry`:
  * Sets the detector geometry to apply acceptance for
  * Options are `4pi` and `LHCb`
  * Default: `LHCb`

* `energy`:
  * Sets the pp collision energy (in TeV) used to get the correct parent 
    kinematics
  * Supported: `7`, `8`, `13` and `14`
  * more types may be added in rootfiles/fonll

* `parent`:
  * Overrides to determined flavour of the parent particle used to get the 
    correct parent kinematics
  * Supported: `b` and `c`
  * More types may be added in [`rootfiles/fonll`][fonll_dir]

* `ptRange`:
  * Define the pT range to generate the parent in.
  * Syntax: `ptRange : <min> <max>`
  * Defaults: 0–300 (for `4pi` geometry) or 0–100 (for `LHCb` geometry)

* `etaRange`:
  * Define the pseudorapidity range to generate the parent in.
  * Syntax: `etaRange : <min> <max>`
  * Defaults: -8–8 for 4pi geometry and 1–6 for LHCb geometry

* `minWidth`:
  * Sets the minimum resonance width (in GeV) to be generated
  * Narrower resonances will be generated with a fixed mass
  * Default: 0.001 GeV

* `maxAttempts`:
  * Sets the maximum number of attempts allowed to generate each decay
  * Default: 1000

* `paramsStable`:
  * Defines the set of parameters to be added to the histograms/tree
    for each stable particle in the decay.
  * Syntax: `paramsStable : <params>`, where `<params>` is a space separated 
    list of the types defined in [the Parameters section][#parameters].
  * Example: `paramsStable : P PT`

* `paramsDecaying`:
  * Defines the set of parameters to be added to the histograms/tree
    for each decaying particle in the decay.
  * Syntax is `paramsDecaying : <params>` where `<params>` is a space separated 
    list of the types defined in [the Parameters section][#parameters].
  * Example: `paramsDecaying : M P PT`

* `paramsTwoBody`:
  * Defines the set of parameters to be added to the histograms/tree
    for each two-body combination of particles in a 3+-body decay.
  * Syntax is `paramsTwoBody : <params>` where `<params>` is a space separated 
    list of the types defined in [the Parameters section][#parameters].
  * Example: `paramsTwoBody : M M2`

* `paramsThreeBody`:
  * Defines the set of parameters to be added to the histograms/tree
    for each three-body combination of particles in a 4+-body decay.
  * Syntax is `paramsThreeBody : <params>` where `<params>` is a space 
    separated list of the types defined in [the Parameters 
    section][#parameters].
  * Example: `paramsThreeBody : M M2`

* `param`:
  * Defines a new parameter to be added to the histograms/tree
  * Syntax is `param : <name> <type> <particles> [TRUE]`, where:
    * `<type>` must be one of the types defined in [the Parameters 
      section][#parameters],
    * `<particles>` is a space separated list of particle indices
    * The optional `TRUE` argument means the parameter will be calculated using 
      the true unsmeared momenta of the particles
  * Example: `param : mSq12 M2 1 2 TRUE`

* `cut`:
  * Applies a cut on a given parameter
  * Syntax is `cut : <param> <type> <min/max> [<max>]`, where
    * `<param>` is the name of a parameter (must be defined using `param`)
    * `<type>` is one of "min", "max", "range" or "veto"
    * `<min>` and/or `<max>` define(s) the cut value(s)

* `shape`:
  * Sets a 1D or 2D PDF to generate events according to
  * Syntax is `shape : <file> <hist> <paramX> [<paramY>]`, where:
    * `<file>` is the file containing the histogram,
    * `<hist>` is the name of the histogram,
    * `<paramX>` is the name of the parameter on the X-axis (must be defined 
      using `param`)
  * For 2D histograms `<paramY>` is the name of the parameter on the Y-axis
  * RapidSim supports TH1F, TH1D, TH2F and TH2D histograms,
  * The dimensionality of the histogram will be inferred from the number of 
    parameters given
  * Defaults to phase-space distribution

### Particle settings

* `name`:
  * A user-friendly name for this particle to be used in variable names
  * Default: Automatically generated unique name of particle

* `smear`:
  * The type of momentum smearing to apply to this particle
  * Supported: `LHCbGeneric`, `LHCbElectron`, `AtlasMuon`, or `AtlasHadron`
    * More types may be defined in config/smear
  * Default: `LHCbElectron` (for electrons/positrons), otherwise `LHCbGeneric`

* `invisible`:
  * Whether the particle should be treated as invisible
  * Invisible particles are not included when determining non-truth parameters
  * A corrected mass parameter will be added automatically if any particles are invisible
  * Default: `true` (for neutrinos), otherwise `false`

* `altMass`:
  * Adds alternative mass hypotheses for this particle
  * Syntax is `altMass : <particles>` where `<particles>` is a space separated 
    list of particle type names as listed in [`config/particles.dat`][particles_file]
  * If alternative hypotheses are defined for multiple particles then 
    parameters will be calculated for all combinations of one or two 
    mis-identifications

## Parameters

* `M`: The invariant mass of the combination of the given particles
* `M2`: The squared invariant mass of the combination of the given particles
* `MT`: The transverse mass of the combination of the given particles
* `E`: The energy of the combination of the given particles
* `ET`: The transverse energy of the combination of the given particles
* `P`: The total momentum of the combination of the given particles
* `PX`: The X momentum of the combination of the given particles
* `PY`: The Y momentum of the combination of the given particles
* `PZ`: The Z momentum of the combination of the given particles
* `PT`: The transverse momentum of the combination of the given particles
* `eta`: The pseudorapidity of the combination
* `phi`: The azimuthal angle of the combination
* `y`: The rapidity of the combination
* `gamma`: The relativistic gamma factor of the combination
* `beta`: The velocity of the combination
* `theta`: The angle between the first two particles in the rest frame of 
  the combination of the remaining particles
  e.g. 1 2 1 3 would give the angle between 1 and 2 in the rest frame 
  of a resonance in `m_13` (i.e. `theta_13` in a 3-body decay)
* `costheta`: The cosine of theta
* `Mcorr`: The corrected mass of the combination of the given particles
  correcting for any invisible particles

## TODO

* Add more particles to [`particles.dat`][particles_file]
* Add derived class for particle specific parameters to simplify RapidParam interface
* Combine `.decay` and `.config` files to simplify interface

## Authors

 * Greig Cowan 2015
 * Dan Craik 2016
 * Matt Needham 2015

[tgenphasespace]: https://root.cern.ch/doc/master/classTGenPhaseSpace.html
[fonll]: http://cacciari.web.cern.ch/cacciari/fonll/
[fonll_dir]: rootfiles/fonll
[particles_file]: config/particles.dat
