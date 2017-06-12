Files in this folder are utilities needed for the Laplacian-DTI Comparison analysis.
Two were written or edited to get around using the Statistics Toolbox and two are 
Matlab File Exchange files that are needed as utilities to run other scripts and 
functions in the Laplacian-DTI Comparison project.



gh_rangesearch 
	-- a workaround to use rangesearch without the Statistics Toolbox
	It uses John D'Errico's ipdm function and manipulates the answer to return the same format as rangesearch.
	Should be able to use rangesearch documentation



ipdm
	-- interpoint distance matrix
	quickly computes a matrix of distances between all points in an array
	written by John D'Errico




YiCao_knnsearch
	-- a version of knnsearch that can be used without the Statistics Toolbox
	version written by Yi Cao and edited by Geoff Handsfield to be consistent with knnsearch




inpolyhedron
	-- tests if points are inside a 3D triangulated surface
	written by Sven Holcombe