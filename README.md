# Welcome to my Exploration of various Computational Geometry Concepts!

### Will include demos of different geometrical algorithms in 2d and 3d space.

---

## *Usage Instructions below Examples*

# RRT-Multi-Joint-Path-Planning

##  Examples
### Single Joint!


https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/429f9e4f-3f77-485d-a020-4bdff71a4e19


## Multi-Joint!


https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/6a87f024-2b1b-4c44-a87e-36415cbf26a4


## My Personal Favorite!


https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/bee5a19b-fdd9-4abb-89ea-5b3ce07e2bc5


## Main Flow
### Begin by drawing a polygonal obstacle, counter-clockwise
### Press n to save current polygon, and begin drawing next
### Once finished drawing obstacles, press n and then s to start drawing robot
### press s to draw a new section of the robot
### press e to choose the goal location
### press r to run RRT
### press p to show path
### press g to go

## q to quit!


---
## OTHER

### t: hide/show tree structure
### c: show all possible conditions

---
## Parameters (Set in CPP file)

### epsilon: controls distance between nodes
### entropyThreshhold: controls the amount the robot can rotate between nodes
### MovementSpeed: speed to animate robot
### SizeRRT: upper limit for amount of nodes in RRT tree

---
---

# Art-Gallery-Guarding

## Examples!

https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/1595e997-aa6a-44ce-a8db-bd5d9e32d6a2

## Usage 

1. Run make command to generate executable
2. Run executable
3. Draw boundaries of shape with left-mouse clicks
4. press 's' to then enter the guard drawing mode
5. click anywhere to add guards
6. press 'r' to let guards move
7. toggle 'v' to show/unshow visibility zones for each of the guards
8. 'q' to quit!

---
---

# Find-Visibility-Graph

## Examples!

https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/cdfc2d31-156c-4742-b59d-f58381833901

## Usage

1. Run make command to generate executable
2. Run executable
3. Draw boundaries of shape with left-mouse clicks
4. Press 'n' to add new boundary (can create as many boundaries as desired)
5. press 's' to then enter the start drawing mode
6. After adding the start, press 'e' to add the end node
7. Press 'r' to run!
8. press 'q' to quit

---
---

# 3D-Hull-Generation

## Examples!

https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/d0e1858b-4782-4b78-b342-6eee483ca243


<img width="428" alt="Screenshot 2023-04-15 at 7 46 48 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/8263bd30-a994-458d-ac9f-97fbdc4f8f15">

<img width="428" alt="Screenshot 2023-04-16 at 12 40 40 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/46b9312a-80f6-403c-a584-e85205783f79">

<img width="450" alt="Screenshot 2023-04-17 at 4 47 09 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/5387d712-3def-45dc-85cd-025c91fb9ae2">

## Usage

1. run "make" command to build project
2. execute the  excecutable with an integer value corresponding to the number of points to generate a hull around!
3. Use x, y, and z to rotate around those respective axes
4. Use q to quit
5. use a to animate the generation process
6. use t to turn on/off geometry around the hull

---
---

# KD-Tree-Mondrian-Generation

## Examples!

<img width="924" alt="Screenshot 2023-02-28 at 3 16 41 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/6239b255-be53-4f02-90be-1527d1d11f2f">

<img width="1001" alt="Screenshot 2023-02-28 at 3 15 16 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/7514a83c-df50-4791-a9c7-0049cb978924">

<img width="915" alt="Screenshot 2023-02-28 at 3 12 42 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/657fd103-55f6-4588-8a8f-331c4d8bfc24">


<img width="1001" alt="Screenshot 2023-02-28 at 3 17 30 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/120df7c4-e069-46bc-92da-0e3f2e44ae0f">


## Usage
1. run "make" command to build project
2. execute the  excecutable with an integer value corresponding to the number of sections to divide the mondrian-style painting into!


---
---

# 2D Hull Generation

## On a given set of points, will generate a 2d hull around them, following the Graham's scan algorithm

## Examples:

<img width="495" alt="Screenshot 2023-09-13 at 1 03 53 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/20caa8a4-8f71-468b-a366-00de5aaa9b83">

<img width="495" alt="Screenshot 2023-09-13 at 1 04 00 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/25c42daf-53ff-487b-b612-7fcd275563bb">

<img width="495" alt="Screenshot 2023-09-13 at 1 04 25 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/aff22197-6387-4ef9-a288-63c65b794b57">

<img width="495" alt="Screenshot 2023-09-13 at 1 04 28 PM" src="https://github.com/DJRGVC/Exploring-Computational-Geometry/assets/42795164/ccae529b-6696-4dea-9ac3-3537e1c40902">

## Usage

1. run "make" command to build project
2. execute the viewPoints excecutable with an integer value corresponding to the number of points to generate
3.  press i to cycle through different shapes, and be amazed as 2d hulls are generated!

---

# Find-Closest-Pair 

## To find the closest pair, we used the following approach:
1. First split grid with k divisions, where k is either the number inputted by the user for number of grid divisions, or is the square root of n, as this should be the optimal number of grid divisions given that the points are randomly distributed, as in this case the number of grid boxes is equal to the number of points.
2. Now, after having decided on the number of grid divisions, we inserted each of the randomly generated points into a 2d array of vectors, where each vector represents the points in one of the k^2 grid boxes. 
3. After having inserted the points into the 2d array, we now call our gridding function to determine which pair of points is the closest.
4. Then, for each point in p, we check the distances from that point to all other points in its grid box. After doing so, we check the distance from that point to other points in adjacent boxes that are within the current minimum distance away from the point, where the current minimum distance is the global minimum distance thus far between two points.
5. If there are no points in a given grid box, we simply skip to the next one.
6. As an added optimization, our gridding approach terminated after a pair of points are found that are a distance one away from eachother, as this is the minimum distance possible in the graph.

## To chose the optimal grid size, we did the following:
* Given n points, we determined that the optimal grid division (value for k) on average was the square root of n, as this would give n total grid boxes, so the relationship between the number of points and the number of grid boxes is one-to-one, which gives the fastest runtime on average.

## Table of our algorithm's runtimes:

|    n	 |	 Naive  	|	 Gridding  	|	
|--------|------------|-------------|	
| 5	  	 |	0.000012	|	0.000004    |
| 10	   |	0.000016	|	0.000005    |
| 50	   |	0.000088	|	0.000015    |
| 100  	 |	0.000320	|	0.000021    |
| 500	   |	0.007676	|	0.000075    |
| 1000	 |	0.030955	|	0.000193    |
| 5000	 |	0.693969	|	0.000043    |
| 10000	 |	2.643573	|	0.000006    |
| 50000	 |	65.636965	|	0.000008    |   
| 100000 |  143.382738| 0.000010    |


