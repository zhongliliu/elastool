from elate import ELATE

def visualize_elastic_anisotropy(
    elastic_tensor,     
    density,
    plot=None,
    elastic_calc=None,
    npoints=100,
    show=True,
    
    #keywords for plot_3D_slice
    normal=(1, 0, 0),
    origin=(0, 0, 0)
):
    """
    This method visualizes the elastic properties
    of a material using a provided elastic tensor.
    """

    rowsList = []

    row = elastic_tensor.shape[0]
    col = elastic_tensor.shape[1]

    for i in range(row):
        columnsList = []
        for j in range(col):
            columnsList.append(round(elastic_tensor[i, j], 3))
        rowsList.append(columnsList)

    elastic_tensor = ELATE(rowsList, density) 

    if plot == "2D":
       # elastic_tensor.print_properties_2D()
        fig = elastic_tensor.plot_2D(
            elastic_calc=elastic_calc, npoints=npoints, apply_to_plot=None, show=show)
    elif plot == "3D":
        #elastic_tensor.print_properties()
        fig = elastic_tensor.plot_3D(
            elastic_calc=elastic_calc, npoints=npoints, show=show)
    elif plot == "3D_slice":
        #elastic_tensor.print_properties()
        fig = elastic_tensor.plot_3D_slice(
            elastic_calc=elastic_calc, npoints=npoints, normal=normal, show=show)

    #elastic_tensor.print_properties()
    #if plot == "2D":
    return fig
    #elif plot == "3D":
    #    return meshes
   # else:
   #     return None


def print_elastic_tensor(elastic_tensor, dim, density):
    """
    This function prints the elastic tensor in a formatted manner.
    """

    rowsList = []

    try:
        row = elastic_tensor.shape[0]
        col = elastic_tensor.shape[1]

        for i in range(row):
            columnsList = []
            for j in range(col):
                columnsList.append(round(elastic_tensor[i, j], 3))
            rowsList.append(columnsList)

        elastic_tensor = ELATE(rowsList, density) 

        if dim == "2D":
            elastic_tensor.print_properties_2D()
        elif dim == "3D":
            elastic_tensor.print_properties()
        else:
            print("Invalid dimension. Please specify '2D' or '3D'.")

    except Exception as e:
        print(f"An error occurred: {e}")


