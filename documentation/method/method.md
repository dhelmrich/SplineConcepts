# Spline Mesh Definition

Spline Meshes: Descriptions of a mesh using splines that can be parametrized.

The basic idea is: using spline definitions, a template is created that can be used to generate a mesh for a plant organ section, such as a meristem-leaf junction. The template is then fitted onto the CPlantBox output geometry, using the mesh parameters as input.

A mesh parameter set $P$ is a set of attributes that define the shape of the mesh. All other properties, aside from these parameters, should be within the template. An example of this would be:

The meristem parameter set $P$ is defined by the following attributes:
- $P_{\text{length}}$: The length of the meristem.
- $P_{\text{radius}}$: The radius of the meristem.
- $P_{\text{angle}}$: The angle of the meristem.
- $P_{\text{blade}}$: The width of the leaf blade.

This parameter set is then used to fit the splines, which have pre-defined points at which they are stretched by parameters or stretched by conditional dependency between each other.
