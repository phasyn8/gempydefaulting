import pyvista as pv
import numpy as np
#import gemgis as gg
#import gempy as gp
import pandas as pd
import re

def surfaces_df(geo_model, color_backgroud_cells=True):
    """
    Output a pandas dataframe with the structural relation data
    
    Inputs: 
    geo_model: (Computed Gempy geo_model)
    color_background_cell flag: (bool) Default is True...

    Outputs pandas DataFrame
    """
    
    df = pd.DataFrame()

    group = geo_model.structural_frame.structural_groups # extract the group and element data
    
    # build empties and key lists 
    row_dict = {}
    keys = ["Group", "Element", "Structural Relation", "Color", "is_active"]
    values = []
    elements = []

    # parse through the structural groups and organize into dictionaries then convert to DataFrames
    for i in range(len(group)):
        elements = group[i].elements
        for j in range(len(elements)):
            values = []
            
            values.append(group[i].name)  # group name
            values.append(elements[j].name) # element name 
            values.append(group[i].structural_relation.name) # structural relation type e.g. (fault, onlap, erode...)
            values.append(elements[j].color) # surface color
            values.append(True) if elements[j].is_active else values.append(False) # is active flag
            for index in range(len(keys)):
                row_dict[keys[index]] = values[index] #build dictionary Row
            df_dict = pd.DataFrame([row_dict]) #convert to dataframe row
            df = pd.concat([df, df_dict], ignore_index=True)
            
    """ This only operates if Background color flag is True"""
    if color_backgroud_cells:
    # Function to check if a value is a valid RGB hex color code
        def is_valid_rgb(value):
            if isinstance(value, str) and re.fullmatch(r"#([A-Fa-f0-9]{6})", value):
                return True
            return False

    # Function to apply background color if the cell value is a valid RGB hex color code
        def color_cells(val):
            if is_valid_rgb(val):
                return f'background-color: {val}'
            return ''

    # Apply the color function to the DataFrame
        styled_df = df.style.applymap(color_cells)

    # Return the styled DataFrame
        return styled_df
    
    # else the uncolored one
    else:
        return df
    

def return_surface_from_geomodel(geo_model, surface=None):

    """
    Gather vertices and faces to create polydata set
    input:  gempy - (geo_model) computed gempy geomodel for extraction
            surface - (int) address of surface to mesh

    
    this is a temporary functions as this is working in gemgis v1.2  "create_depth_map_from gempy" or something like that
    
    """
# collect vertices
    vertices = geo_model.input_transform.apply_inverse(geo_model.solutions.raw_arrays.vertices[surface])
# collect faces
    faces = np.hstack(np.pad(geo_model.solutions.raw_arrays.edges[surface], ((0, 0), (1, 0)), 'constant', constant_values=3))
 
# combine into Polydata Mesh
    mesh = pv.PolyData(vertices, faces)

# Build scaler for the mesh
    mesh['Depth [m]'] = mesh.points[:, 2]

    return mesh



def gempydefaulting(geo_model, clip_scalar, clip): #ver 003
    '''
    Method removes all the block stepping artifacts with faults in gempy models
    Output is a Pyvista Plotter with preserved coloration from the input gempy model.
    Methodology, divide the input model into lithological surfaces and fault surfaces,
    clip each surface by each fault with a distance defined by the clip_scalar, this
    should roughly the diagonal dimention of your gempy model. Each edge is then 
    isolated and extruded to the fault and again clipped with the fault to create a 
    clean interface between each surface and each fault.
    
    **** NOW WORKING with GEMPY Ver 3.0 ****
    '''
    #modeldf = gp.get_data(model, 'surfaces').df

    #Isoloate fault surfaces
    #faultsdf = modeldf.where(modeldf['isFault']==True)
    #faultsdf = faultsdf.dropna()

    #elements = model.structural_frame.elements_names

    #Isolate lithological surfaces 
    #surfacesdf = modeldf.where(modeldf['isFault']==False)
    #surfacesdf = surfacesdf.dropna()
    
    #Build list of Series Surfaces and faults
    surfaces = surfaces_df(geo_model, color_backgroud_cells=False)


    surfaces_list = surfaces[surfaces["Structural Relation"]!="FAULT"].index.tolist()
    #model.surfaces.df[model.surfaces.df['isFault']==False]['surface'].to_list()[:-1]
    
    faults_list = surfaces[surfaces["Structural Relation"]=="FAULT"].index.tolist()
    #model.surfaces.df[model.surfaces.df['isFault']==True]['surface'].to_list()
    #cali_df
    #faults_list, surfaces_list

    dict_surfaces = {}
    fault_surfaces = {}

    for i, idx in zip(range(len(surfaces_list)),surfaces_list):
        dict_surfaces[i] = return_surface_from_geomodel(geo_model, surface=idx)


    for i, idx in zip(range(len(faults_list)),faults_list):
        fault_surfaces[i] = return_surface_from_geomodel(geo_model, surface=idx)

    print("faults " + str(faults_list) + "   surfaces  " + str(surfaces_list))
    p = pv.Plotter(notebook=True)
    
    for surface in dict_surfaces:
            layer = surface
            print('**********' + str(layer) + '************')
            working_surface = dict_surfaces[layer]
            working_surface.compute_normals(cell_normals=False, inplace=True)
            flipZ = np.array((1,1,-1))
            for surface in fault_surfaces:

                fault = surface

                print('**********' + str(fault) + '************')
                #add normals to fault surface object
                faults = fault_surfaces[fault]#.delaunay_2d()
                # array of normals
                faults.compute_normals(cell_normals=True, flip_normals=False, inplace=True)
                fault_norm = np.average(faults['Normals'], axis=0)

                #faults.compute_normals(cell_normals=True, flip_normals=False, inplace=True)
                #fault_norm_neg = np.average(faults['Normals'], axis=0)
                # Average vector of all fault surface normals 
                print('fault Norm ' + str(np.around(fault_norm, decimals=2)))
                #print('fault norm neg ' + str(np.around(fault_norm_neg, decimals=2)))

                if fault_norm[0] >= 0:
                    print('Fault normal >= 0')
                        #Clipping each surface with clip_scalar distance from faults
                    left_meshes = working_surface.clip_surface(fault_surfaces[fault], invert=False, value=1*clip_scalar)

                            #finding left edges

                    left_edges = left_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault], invert=True, value=1*(clip_scalar+1))
                    #left_edges.compute_normals(cell_normals=True, inplace=True)

                    left_edges_norm = np.average(left_edges['Normals'], axis=0)*flipZ
                            #average left_edges normals

                    print('left_edges '+' average normal ' + str(np.around(left_edges_norm, decimals=2)))

                            #computing vector for extrusion, PV Extrusion method only accepts one extrusion vector so
                            #squared vector enhances directionality...
                    dirL = (-3*clip_scalar*(left_edges_norm + fault_norm)**2)

                    print('Left ' + str(layer) + ' scalar dir ='+ str(np.around(dirL, decimals=2)))

                            #Extruding just the edge to extend it past the fault  
                    left_meshes_ext = left_edges.extrude(dirL).clip_surface(fault_surfaces[fault], invert=False)

                    right_meshes = working_surface.clip_surface(fault_surfaces[fault], invert=True, value=-1*clip_scalar)

                    right_edges = right_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault], invert=False, value=-1*(clip_scalar+1))
                    right_edges.compute_normals(cell_normals=True, inplace=True)
                    right_edges_norm = np.average(right_edges['Normals'], axis=0)*flipZ
                        #right_edges

                    print('Right edges avg Norm ' + str(np.around(right_edges_norm, decimals=2)))
                    dirR = (3*clip_scalar*(right_edges_norm + fault_norm)**2)
                    print('Right ' + str(layer) + ' scalar dir ='+ str(np.around(dirR, decimals=2)))
                    right_meshes_ext = right_edges.extrude(dirR).clip_surface(fault_surfaces[fault], invert=True)
                        #left_meshes_ext = left_edges.extrude(left_edges['Normals']).clip_surface(fault_surfaces['Main_Fault'][0], invert=True)

                    p.add_mesh(right_meshes_ext, color=surfaces.iloc[layer]["Color"])
                            #Adding each mesh and extension to the Plotter


                    p.add_mesh(left_meshes_ext, color=surfaces.iloc[layer]["Color"])
                    working_surface = left_meshes + right_meshes
                    #working_surface += right_meshes 
                    #working_surface += left_meshes
                    #p.add_mesh
                elif fault_norm[0] < 0:
                    print('Fault normal < 0')
                    #fault_norm = faults.compute_normals(cell_normals=True, flip_normals=True)

                    # Average vector of all fault surface normals 

                    #fault_norm = np.average(fault_norm['Normals'], axis=0)
                    left_meshes = working_surface.clip_surface(fault_surfaces[fault], invert=True, value=-1*clip_scalar)
                    #left_meshes.compute_normals(cell_normals=False, inplace=True)
                            #finding left edges
                    left_edges = left_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault], invert=False, value=-1*(clip_scalar+1))
                    if np.isnan(left_edges['Normals']) or None: 
                        pass 
                    
                    else:

                        print('left edges ' + str(left_edges['Normals']))


                        left_edges_norm = np.average(left_edges['Normals'], axis=0)
                                #average left_edges normals

                        print('left_edges '+' average normal ' + str(np.around(left_edges_norm, decimals=2)))
                                #computing vector for extrusion, PV Extrusion method only accepts one extrusion vector so
                                #squared vector enhances directionality...
                        directionL = (-3*clip_scalar*(left_edges_norm - fault_norm)**2)
                        print('Left ' + str(layer) + ' scalar direction ='+ str(np.around(directionL, decimals=2)))
                                #Extruding just the edge to extend it past the fault  
                        left_meshes_ext = left_edges.extrude(directionL).clip_surface(fault_surfaces[fault], invert=True)



                        right_meshes = working_surface.clip_surface(fault_surfaces[fault], invert=False, value=1*clip_scalar)
                        #right_meshes.compute_normals(cell_normals=False, inplace=True)

                        right_edges = right_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault], invert=True, value=1*(clip_scalar+1))

                        #print(right_edgeright_edges = right_meshes.compute_normals(cell_normals=False)s['Normals'])
                        right_edges_norm = np.average(right_edges['Normals'], axis=0)
                            #right_edges

                        print('right edges avg normal' + str(np.around(right_edges_norm, decimals=2)))
                        directionR = (3*clip_scalar*(right_edges_norm - fault_norm)**2)
                        print('Right ' + str(layer) + ' scalar direction ='+ str(np.around(directionR, decimals=2)))
                        right_meshes_ext = right_edges.extrude(directionR).clip_surface(fault_surfaces[fault], invert=False)
                            #left_meshes_ext = left_edges.extrude(left_edges['Normals']).clip_surface(fault_surfaces['Main_Fault'][0], invert=True)
                        working_surface = left_meshes + right_meshes
                        p.add_mesh(right_meshes_ext, color=surfaces.iloc[layer]["Color"])
                                #Adding each mesh and extension to the Plotter


                        p.add_mesh(left_meshes_ext, color=surfaces.iloc[layer]["Color"])
                        #print('joning meshes')
                    #working_surface += left_meshes
                    #working_surface += right_meshes
            #p.add_mesh(right_meshes, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])          
            #p.add_mesh(left_meshes, color='blue'.to_list()[0])
            p.add_mesh(working_surface, color=surfaces.iloc[layer]["Color"])
            for surface in fault_surfaces:
                p.add_mesh(fault_surfaces[surface], color=surfaces.iloc[surface]["Color"])

    p.show_grid(color='black')
    p.set_background(color='white')
    p.set_scale(1,1,1)

    p.show()
