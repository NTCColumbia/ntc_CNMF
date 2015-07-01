# ca_source_extraction

The code implements a method for simultaneous source extraction and spike inference from large scale calcium imaging movies. The code is suitable for the analysis of somatic imaging data. Implementation for the analysis of dendritic/axonal imaging data will be added in the future. Depending of the chosen method the algorithm may need the constrained deconvolution method maintained separately in https://github.com/epnev/constrained-foopsi 

The algorithm is presented in more detail in

Pnevmatikakis, E. A., Gao, Y., Soudry, D., Pfau, D., Lacefield, C., Poskanzer, K., ... & Paninski, L. (2014). A structured matrix factorization framework for large scale calcium imaging data analysis. arXiv preprint arXiv:1409.2903. http://arxiv.org/abs/1409.2903

function name                           | description 
----------------------------------------|-----------------------------------
demo_script.m                           | wrapper code <br />
update_spatial_components.m             | update spatial components given temporal components and data <br />
update_temporal_components.m            | update temporal components given spatial components and data <br />
merge_ROIs.m                            | merge spatially overlapping components that are temporally correlated <br />
utilities/arpfit.m                      | estimation of noise level for every pixel and global time constants <br />
utilities/com.m:                        | calculation of the center of mass of each component <br />
utilities/correlation_image.m           | calculates the correlation image of the movie <br />
utilities/graph_connected_comp.m        | finds the connected components in a graph <br />
utilities/greedyROI2d.m                 | Greedy method for initializing the spatial and temporal components <br />
utilities/interp_missing_data.m         | Filling in missing data using linear interpolation <br />
utilities/lars_regression_noise.m       | solve a basis pursuit denoising problem using the LARS algorithm <br />
utilities/make_G_matrix.m               | construct a convolution/deconvolution matrix <br />
utilities/make_patch_video.m            | construct a video that displays the results of the algorithm <br />
utilities/order_ROIs.m                  | order found components based on their temporal activation and spatial size <br />
utilities/plain_foopsi.m                | projection of fluorescence onto the cone formed by the indicator dynamics 
utilities/plot_contours.m               | contour plot of found components and creation of a json file <br />
utilities/tiff_reader.m                 | loading a tiff stack into matlab <br />
utilities/view_patches.m                | plotting of each found component and its temporal activation <br />


License
=======

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
