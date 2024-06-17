# inviscid-fluid-simulation

This project uses the finite volume method to solve the Euler equations used to describe the time evolution of an invicid fluid. An animation is produced displaying the Kelvin Helmholtz intability patterns that develop due to the chosen initial conditions. These initial conditions describe two fluids of differing densities running in opposite directions with small perturbations at their interfaces.

# Usage

kelvin_helmholtz.cpp is used to produce output density data that is then used by highres-fluid-animation.ipynb to create images and animations.

All of these files used should be fairly easy using the command line and a Jupyter Notebook. Be sure to allow the files to access eachother by placing the .cpp files in the same directory. It will also be necessary to provide appropriate file paths in highres-fluid-animation.ipyb.

# Animations

HD_Fluid_Anim.mp4 provides an example of the output of this project's pipeline.
