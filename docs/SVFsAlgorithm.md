# Single Vertex Fitting on straight tracks (SVFs)

> :bulb: **tl;dr**: <br>
This fitter is a single secondary vertex fitter. It is a straight line fitter that iteratively minimizes a $\chi^2$ defined by the sum of all vertex-to-track distances.

The SVFs (**S**ingle **V**ertex **F**inder on **s**traight tracks) algorithm implemented in this code (`modules.singleVertexFitter`) is described in details in [this book](https://link.springer.com/book/10.1007/978-3-030-65771-0), in paragraph `8.1.1.1`, and the notation is kept as close as possible in the code. The fitter builds the $\chi^2$ as the sum over all the distances between the fitted vertex and the tracks, assuming tracks to be straight lines. Then, it minimizes the $\chi^2$ iteratively using the Newton-Raphson method, until stability (meaning that the the vertex update is within a certain tolerance) or maximum iteration number is reached.
 
## Iteration Steps
Given the tracks origins $r_i$ and versors $a_i$, and defining $t$ as the iteration index and $i$ as the tracks index:
- Build the jacobian of the least square loss $D$ for each track origin $r_i$ as:
$$\frac{\partial D_{i,t}}{\partial r_i} = 2(a_{i,1} \eta_{i,10} - a_{i,2} \eta_{i,02}, a_{i,2} \eta_{i,21} - a_{i,0} \eta_{i,10},a_{i,0} \eta_{i,02} - a_{i,1} \eta_{i,21})^\top$$

- Build the jacobian of the least square loss $D$ for each track versor $a_i$ as:
$$\frac{\partial D_{i,t}}{\partial a_i} = 2(d_{i,2} \eta_{i,02} - d_{i,1} \eta_{i,10}, d_{i,0} \eta_{i,10} - d_{i,2} \eta_{i,21},d_{i,1} \eta_{i,21} - d_{i,0} \eta_{i,02})^\top$$

with the auxiliary variables:
$$d_{i,k } = r_{i,k} - v_{t,k}$$ 
where $v$ is the fitted vertex at iteration $t$;
$$\eta_{i,jk} = a_{i,j}d_{i,k} - a_{i,k}d_{i,j}$$

- Build the Jacobian of the track as:
$$J_i = (\frac{\partial D_{i,t}}{\partial r_i} ,\frac{\partial D_{i,t}}{\partial a_i})^\top$$

- Evaluate the tracks' weight as:
$$\sigma^2_i = J^\top V_i J$$
where $V$ is the covaraince matrix of the track. 
> **Note**
To control stability a small quantity is added to $\sigma^2$.

- Evaluate the total gradient of the square loss as:
$$\nabla D = - \sum\limits_{i} \frac{1}{\sigma^2_i} \frac{\partial D_{i,t}}{\partial r_i}$$

- Build the hessian of each track as:
$$diag(H_i) = (a_{i,1}^2 + a_{i,2}^2,a_{i,0}^2 + a_{i,2}^2,a_{i,0}^2 + a_{i,1}^2)^\top$$
$$H_{i,01} = H_{i,10} = -a_{i,0}a_{i,1}$$
$$H_{i,02} = H_{i,20} = -a_{i,0}a_{i,2}$$
$$H_{i,12} = H_{i,21} = -a_{i,1}a_{i,2}$$


- Calculate the laplacian of the square loss as:
$$\nabla^2 D = \sum\limits_{i} \frac{2}{\sigma^2_i} H_i$$

- Update the vertex as:
$$v_{t+1} = v_t - (\nabla^2 D )^{-1} \nabla D$$

## Iterative minimization

The $\chi^2$ of the fit is minimized iteratively with the Newton-Raphson method. The iterative minimization ends when the maximum iteration number is reached or when the vertex update is smaller than a tolerance threshold (`eps` in the class constructor):
$$|v_{t+1} - v_t|<\varepsilon$$

> **Note**
To avoid stopping the fit if one iteration is randomly stuck, the fit is actually stopped when the stability criteria is hit 5 times in a row.