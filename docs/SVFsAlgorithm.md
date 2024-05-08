# Single Vertex Fitting on straight tracks (SVFs)

> :bulb: **tl;dr** <br>
This fitter is a single secondary vertex fitter. It is a straight line fitter that iteratively minimizes a $\chi^2$ defined by the sum of all vertex-to-track distances.

The SVFs (**S**ingle **V**ertex **F**inder on **s**traight tracks) algorithm implemented in this code (`modules.singleVertexFitter`) is described in details in [this book](https://link.springer.com/book/10.1007/978-3-030-65771-0), in paragraph `8.1.1.1`, and the notation is kept as close as possible in the code. The fitter builds the $\chi^2$ as the sum over all the distances between the fitted vertex and the tracks, assuming tracks to be straight lines. Then, it minimizes the $\chi^2$ iteratively using the Newton-Raphson method, until stability (meaning that the the vertex update is within a certain tolerance) or maximum iteration number is reached.

## $\chi^2$ definition 
Given a track, represented as straight line by its origin $\boldsymbol{r}_i$ and its versor $\boldsymbol{a}_i$, and a point in space $\boldsymbol{v}$, their squared distance is equal to:
$$D^2_i(\boldsymbol{v}) = [(\boldsymbol{r}_i - \boldsymbol{v})\wedge \boldsymbol{a}_i]^2$$
The $\chi^2$ of the fit is obtained by summing the distances between the vertex and each track, and weighting by the track's variance:
$$\chi^2(\boldsymbol{v}) = \sum\limits_{i=1}^{N_{tracks}} \frac{D^2_i(\boldsymbol{v})}{\sigma_i^2}$$
This $\chi^2$ is iteratively minimized to perform the fit. In the following we will call the $\chi^2$ as $\mathcal{S}(\boldsymbol{v})$, to keep the notation coherent with the book.

## Iteration Steps
Given the tracks origins $\boldsymbol{r}_i$ and versors $\boldsymbol{a}_i$, and defining $t$ as the iteration index and $i$ as the tracks index:
- Build the partial derivative of the least square loss $D^2_i$ for each track origin $r_i$ as:
$$\frac{\partial D^2_{i,t}}{\partial \boldsymbol{r}_i} = 2(a_{i,1} \eta_{i,10} - a_{i,2} \eta_{i,02}, ~a_{i,2} \eta_{i,21} - a_{i,0} \eta_{i,10},~a_{i,0} \eta_{i,02} - a_{i,1} \eta_{i,21})^\top$$

- Build the partial derivative of the least square loss $D^2_i$ for each track versor $\boldsymbol{a}_i$ as:
$$\frac{\partial D^2_{i,t}}{\partial \boldsymbol{a}_i} = 2(d_{i,2} \eta_{i,02} - d_{i,1} \eta_{i,10}, ~d_{i,0} \eta_{i,10} - d_{i,2} \eta_{i,21},~d_{i,1} \eta_{i,21} - d_{i,0} \eta_{i,02})^\top$$

with the auxiliary variables:
$$d_{i,k } \triangleq r_{i,k} - v_{t,k}$$ 
where $v$ is the fitted vertex at iteration $t$;
$$\eta_{i,jk} \triangleq a_{i,j}d_{i,k} - a_{i,k}d_{i,j}$$

- Build the Jacobian of the track as:
$$\boldsymbol{J}_i = (\frac{\partial D^2_{i,t}}{\partial \boldsymbol{r}_i} ,\frac{\partial D^2_{i,t}}{\partial \boldsymbol{a}_i})^\top$$

- Evaluate the tracks' weight (variance) as:
$$\sigma^2_i = \boldsymbol{J}^\top \boldsymbol{V}_i \boldsymbol{J}$$
where $\boldsymbol{V}$ is the covaraince matrix of the track's parameters. 
> :memo: **Implementation detail** <br>
To control stability a small quantity is added to $\sigma^2$. Furthermore, only the diagonal of the covariance matrix is available in the `H5` file used in the fit.

- Evaluate the total gradient of the $\chi^2$ as:
$$\nabla \mathcal{S} = - \sum\limits_{i} \frac{1}{\sigma^2_i} \frac{\partial D^2_{i,t}}{\partial \boldsymbol{r}_i}$$

- Build the hessian of each track as:
$$diag(\boldsymbol{H}_i) = (a_{i,1}^2 + a_{i,2}^2,a_{i,0}^2 + a_{i,2}^2,a_{i,0}^2 + a_{i,1}^2)^\top$$
$$H_{i,01} = H_{i,10} = -a_{i,0}a_{i,1}$$
$$H_{i,02} = H_{i,20} = -a_{i,0}a_{i,2}$$
$$H_{i,12} = H_{i,21} = -a_{i,1}a_{i,2}$$


- Calculate the hessian of the $\chi^2$ as:
$$\nabla^2 \mathcal{S} = \sum\limits_{i} \frac{2}{\sigma^2_i} \boldsymbol{H}_i$$

- Update the vertex as:
$$\boldsymbol{v}_{t+1} = \boldsymbol{v}_t - (\nabla^2 S)^{-1} \nabla S$$

## Iterative minimization

The $\chi^2$ of the fit is minimized iteratively with the Newton-Raphson method. The iterative minimization ends when the maximum iteration number is reached or when the vertex update is smaller than a tolerance threshold (`eps` in the class constructor):
$$||\boldsymbol{v}_{t+1} - \boldsymbol{v}_t||<\varepsilon$$

> :memo: **Implementation detail**<br>
To avoid stopping the fit if one iteration is randomly stuck, the fit is actually stopped when the stability criteria is hit 5 times in a row.

When the fit stops at the iteration $T$, the $\chi^2$ of the fit is obtainable as:
$$\chi^2(\boldsymbol{v}_T) = \sum\limits_{i=1}^{N_{tracks}} \frac{D^2_i(\boldsymbol{v}_T)}{\sigma_i^2}$$
and the covariance matrix of the vertex is obtainable as:
$$\left. \nabla^2 \mathcal{S}\right|_{t=T}$$
