# Closure of Meshes

This is a transcription of my blackboard notes.

## Definition

A mesh $M$ consists of a set of splines $\lbrace S_1, \cdots, S_2\rbrace$ and a set of properties $P$.


## Closure

We only look at a section. If there is something about the mesh that is not closed, but individually they are closed, then the mesh just consists of two meshes.

**Perpentidularity**: $`S_1 \perp S_2$ if $\exists t \in \mathscr{D}_\mathbb{R}(S_1) \wedge \exists t' \in \mathscr{D}_\mathbb{R}(S_2) : S_1(t) = S_2(t')`$.

**Closure**: A mesh is closed if $\forall S \in M : \exists S' \in M : S \perp S'$.

This notion of closure is very loose and not really a classical notion of closure, but the important part is that the connection of two splines (probably explicit, because the system should not be bothered to have to find intersections) is sufficient to produce triangles.

**Extrusion**: If $S \in M$ is being extruded using a base shape $G$, then $S$ can be closed using $G$, by connecting copies of $G$ that are being placed along the spline. This is the geometrization version that produces plant meshes in [CPlantBox](https://github.com/Plant-Root-Soil-Interactions-Modelling/CPlantBox/pull/101) currently.

**Ring Closure**: if $S \in M$ is a ring, then it is closed by definition. This leads to undefined behavior, but visual feedback is given to the user. A ring is closed to a mesh by connecting a linspace between the first and last point of the ring.

## Geometrization

1) Find a global grid that can be used across the section.
2) The geometrization can be done if all splines are closed and not degenerate.
3) Choose a primary spline per segment. This primary spline should have the most intersections.
4) Cut all perpendicular splines at the intersection points with the primary spline.
5) Connect the cut splines with the primary spline using SplineConcepts::Geomtry::fillbetween.
6) Repeat for all splines.
