# p1-closest-try-angle

### Project by Brian Liu and Daniel Grant
    a brief decsription of how you chose the grid size and why, and what’s the expected running time assuming the points are uniformly distributed
    the table with the experimental running times of your two algorithms, for various values of n

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
