Goal
----
Statismo-shaperegistration contains example programs that illustrate how [statismo](http://github.com/statismo/statismo) can be used for shape registration using the
idea of Gaussian process registration, as described in the [paper](http://gravis.cs.unibas.ch/publications/2013/MLMI-Luethi_etal-2013.pdf).

Installation
------------
Buidling and installation is done using the standard CMake workflow.

Usage
-----
This usage example uses the example femur surfaces that are provided in the data directory.

In a first step, we choose a reference shape and build a deformation model using this reference with ```buildgpmodel```:

```
buildgpmodel data/VSD001_femur.vtk 50 50 100 femurmodel.h5
```
This build a Gaussian Process model with the name femurmodel.h5, where a Gaussian Kernel of bandwith 50 and scale 50 is approximated using 100 basis functions
(see the [paper](http://gravis.cs.unibas.ch/publications/2013/MLMI-Luethi_etal-2013.pdf) for details.)
The model and the variations can be visualized using the [statismo modelviewer](http://statismo.cs.unibas.ch/statismo-ui.jar).

In the next step, the model is fitted (registered) with the target shape.

```
shaperegistration femurmodel.h5  data/VSD002_femur.vtk 0.01 fitting-result.vtk projection.vtk

```
This command performs a fitting to the target surface VSD002_femur.vtk using a regularization weight of 0.01 and stores the fitting result
and a projection on the target surface. (Note, that it is assumed that the target surface is already aligned with the model.)

In order to force correspondence, it is also possible to specify a number of landmark points that should be matched
```
shaperegistration femurmodel.h5  data/VSD002_femur.vtk 0.01 fitting-result.vtk projection.vtk VSD001-lm.csv VSD002-lm.csv 1

```
The last parameter here specifies the inaccuracy of the landmark placements (it reflects the variance of a gaussian noise model).


Remark
-------
It is not clear yet if the code that is posted here will be maintained and updated for future statismo version.


Licence
-------
This code is licensed under the terms of the BSD license.


Ackowledgment
-------------
The femur datasets that are part of this example have been provided by the [VSD project](https://www.virtualskeleton.ch/).
