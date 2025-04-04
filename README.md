# LEPL1110 - Finite Element Project: Bridge vs. Tank

<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/made-with-c.svg" alt="Made with C"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/contains-tasty-spaghetti-code.svg" alt="Contains Tasty Spaghetti Code"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/built-with-love.svg" alt="Built With Love"></a>
<a href="https://forthebadge.com"><img src="https://forthebadge.com/images/badges/60-percent-of-the-time-works-every-time.svg" alt="60% of the time, works every time"></a>

## The Saga of Yevgeny's Bridge Problem...

Meet Yevgeny. Yevgeny has a T-80 battle tank and a penchant for, let's say, energetic infrastructure testing. While using a certain beam bridge over the Dnieper for target practice, his expert marksmanship (with a 125mm explosive shell, no less) takes out the central beam. Boom.

To his surprise, the bridge, though wounded, defiantly remains standing. Before Yevgeny can deliver a coup de grâce, his lieutenant radios in: "Forget the fireworks, Yevgeny, we need to cross that river and visit the next village!"

Oops. The critical question arises: Can this spectacularly damaged bridge handle the hefty weight of Yevgeny's T-80? That's what this project aims to find out (through the magic of finite element analysis, of course!).

## Table of Contents

- The Saga of Yevgeny's Bridge Problem...
- Project Overview
- Code Structure: Organized Chaos
- Dependencies: The Necessary Evils
- Getting it Running: Compilation Rituals
- How to (Try to) Cross the Bridge: Usage
- Future Dreams & Known Quirks

## Project Overview

This repository contains our code for the LEPL1110 Finite Element Method project. We simulate the structural integrity of a damaged bridge under the load of a tank. It's built primarily in C, featuring a custom FEA solver and OpenGL for visualization.

## Code Structure: Organized Chaos

We've structured the code into two main parts:

- **The Main Project (Root Directory)**: This is the whole shebang. It handles mesh generation (using gmsh), rendering the bridge and tank (via glfw, glad, and our own lightweight graphics code found in src/gui), loading assets (stb_image.h), and logging useful (or cryptic) messages (log.c). It depends on the solver library below.

- **The Standalone Solver (solver/ directory)**: Nestled inside the solver directory is the heart of the FEA simulation. It's an independent library, crafted from scratch with C, containing our finite-element method implementation. Not any dependencies here!

## Dependencies: The Necessary Evils

We tried to keep things contained:

- **Vendored**: Most minor dependencies (glfw, glad, log.c, stb_image.h) are included directly in the repository for your convenience.
- **gmsh**: The one major external dependency. CMake should automagically download the SDK for you.
  - Murphy's Law Clause: If https://gmsh.info happens to be down (like it is right now while I'm writing those lines at 2 AM on April 4th), don't worry. It's OK.  You can manually download the SDK and extract it to `gmsh/gmsh-4.13.1-Linux64-sdk` at the project root. CMake should then detect it and let the build proceed. (Fun fact ! You can download the SDK from the Wayback Machine !)

## Getting it Running: Compilation Rituals

Follow the sacred CMake incantations:

To build and run the full bridge simulation project:
```bash
mkdir build
cmake -B build
make -C build
./build/bridge
```
⚠️ WARNING: You MUST run the bridge executable from the root directory of the project (the one containing this README). It needs to find shaders and assets via relative paths. Running it from the build directory will likely result in sadness, crashes, or summoning unintended demons. (If it does work from the build directory... well, you're a wizard, Harry.)

To build only the standalone solver library:
```bash
cd solver
mkdir build
cmake -B build
make -C build
./build/solver-bin
```
This will also create a static library solver/build/libsolver.a if you wish to link against it elsewhere.

## How to (Try to) Cross the Bridge: Usage

Once the simulation is running:

- **Tank Movement**: Use the Arrow Keys to carefully drive Yevgeny's tank across the bridge. Pray the bridge holds up!
- **Target Practice**: Use your Mouse to aim and click to simulate firing the cannon. See the stress distributions change in real-time !

Shortcuts:

- We had grand plans for implementing helpful shortcuts.
- TODO: Actually implement some shortcuts. Please, benevolent spirit of project deadlines, grant us more time! (If you're reading this, the spirit was probably busy elsewhere.)

## Future Dreams & Known Quirks

- Clean up the legacy mesh code in src/mesh.
- Implement those elusive shortcuts.
- Maybe add more tanks? Or different types of bridges? (Scope creep, hello!)

Enjoy exploring the precarious situation of Yevgeny and his bridge! We hope you (and the grader!) appreciate the simulation and the character of the code.
