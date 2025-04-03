# LEPL1110 - Projet d'éléments finis

<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/made-with-c.svg" alt="Made with C"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/contains-tasty-spaghetti-code.svg" alt="Contains Tasty Spaghetti Code"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/built-with-love.svg" alt="Built With Love"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg" alt="60% of the time, works every time"></a>

## Un peu de contexte...

Yeyvegeny, soucieux de détruire au maximum l’infrastructure du pays qu’il envahit, utilise un pont à
poutres qui traverse le Dniepr comme cible d’entrainement pour son char de combat T-80. Étant un
artilleur expérimenté, Yevgeny touche la poutre centrale en plein milieu avec un obus 125mm
explosif et, à sa surprise, malgré l’effondrement de la poutre centrale, le pont tient bon. Avant de
pouvoir tirer une seconde fois il reçoit un ordre de son lieutenant : ils doivent traverser le fleuve et
attaquer le village le plus proche. Oups! Le pont, dans son état endommagé, pourra-t-il supporter le
poids du char de Yeyvgeny?

## Table of contents

- [Code Structure](#code-structure)

## Code structure

The code is divided into two parts. The project is the repository as a whole, but inside
the `solver` directory you will find a independant library containing a finite-element method solver implemented from scratch without any dependencies.
On the other-hand, the project depends on `gmsh` to generate meshes, as well as `glfw` and `glad` to render them using OpenGL 3.3 standard. We are also using ultra-lightweight libs such as `log` and `stb_image.h` to respectively..log and load images from png files.

All dependencies are vendored, excepted gmsh. However, CMake should download it automagically for you :D
If, just like us at the time of writing this paragraph (4 of April, 2am) the gmsh website is down for no reason and you are stuck, it's ok. don't worry. it's still doable using good old way. Simply extract gmsh sdk into `gmsh/gmsh-4.13.1-Linux64-sdk` at the root of the project, and cmake should detect it and let you continue grading this project :-) (oh and btw you can download the gmsh sdk from the Wayback Machine, nice to know)

To talk a little bit about the code structure, we decided to re-implement a lightweight graphic library using latest OpenGL specs, and its code can be found inside `src/gui` folder.

All the mesh generation stuff can be found inside `src/mesh`, as well as some legacy code used to parse txt mesh files. TODO: clean this up if we do have time. If you can see this text written here, then we didn't.

## Compilation

To run the project, it's the standard cmake procedure:

    mkdir build
    cmake -B build
    make -C build
    ./build/bridge
    
WARNING: The code **MUST** be run from the root of the project, as it will try to load files such as shaders or assets by relative path. If you try running from build folder, then it should not work. ̣̣~And if it does work, then you're a wizard, Harry.~

However, if you just want to build the solver, then let's meet inside the `solver` directory.
Then it's pretty much the same steps:

    mkdir build
    cmake -B build
    make -C build
    ./build/solver-bin

A static library will also be created, named `build/libsolver.a`.

## Usage

The usage is relatively simple. You can use arrow keys to move the **tank** around the bridge. And with your mouse, you can simulate a canon, shooting at the bridge.

Some shortcuts are also implemented such as:

TODO: make some shortcuts. plz god give us more time.

