import pyvista as pv
import numpy as np
import gemgis as gg

def gempydefaulting(model, clip_scalar, clip): #ver 002
    '''Method removes all the block stepping artifacts with faults in gempy models
    Output is a Pyvista Plotter with preserved coloration from the input gempy model.
    Methodology, divide the input model into lithological surfaces and fault surfaces,
    clip each surface by each fault with a distance defined by the clip_scalar, this
    should roughly the diagonal dimention of your gempy model. Each edge is then 
    isolated and extruded to the fault and again clipped with the fault to create a 
    clean interface between each surface and each fault.
    '''
    modeldf = gp.get_data(model, 'surfaces').df

    #Isoloate fault surfaces
    faultsdf = modeldf.where(modeldf['isFault']==True)
    faultsdf = faultsdf.dropna()

    #Isolate lithological surfaces 
    surfacesdf = modeldf.where(modeldf['isFault']==False)
    surfacesdf = surfacesdf.dropna()
    
    #Build list of Series Surfaces and faults
    series_list = surfacesdf['series'].unique().to_list()

    surfaces_list = model.surfaces.df[model.surfaces.df['isFault']==False]['surface'].to_list()[:-1]
    
    faults_list = model.surfaces.df[model.surfaces.df['isFault']==True]['surface'].to_list()
    
    #type(fault)
    #print(surfaces)
    
    #use gemgis to create dictionary of PolyData surfaces
    dict_surfaces = gg.visualization.create_depth_maps_from_gempy(model, surfaces=surfaces_list)
    fault_surfaces = gg.visualization.create_depth_maps_from_gempy(model, surfaces=faults_list)
    p = pv.Plotter(notebook=False)
    
    for surface in dict_surfaces:
            layer = surface
            print('**********' + str(layer) + '************')
            working_surface = dict_surfaces[layer][0]
            working_surface.compute_normals(cell_normals=False, inplace=True)
            flipZ = np.array((1,1,-1))
            for surface in fault_surfaces:

                fault = surface

                print('**********' + str(fault) + '************')
                #add normals to fault surface object
                faults = fault_surfaces[fault][0]#.delaunay_2d()
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
                    left_meshes = working_surface.clip_surface(fault_surfaces[fault][0], invert=False, value=1*clip_scalar)

                            #finding left edges

                    left_edges = left_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault][0], invert=True, value=1*(clip_scalar+1))
                    #left_edges.compute_normals(cell_normals=True, inplace=True)

                    left_edges_norm = np.average(left_edges['Normals'], axis=0)*flipZ
                            #average left_edges normals

                    print('left_edges '+' average normal ' + str(np.around(left_edges_norm, decimals=2)))

                            #computing vector for extrusion, PV Extrusion method only accepts one extrusion vector so
                            #squared vector enhances directionality...
                    dirL = (-3*clip_scalar*(left_edges_norm + fault_norm)**2)

                    print('Left ' + str(layer) + ' scalar dir ='+ str(np.around(dirL, decimals=2)))

                            #Extruding just the edge to extend it past the fault  
                    left_meshes_ext = left_edges.extrude(dirL).clip_surface(fault_surfaces[fault][0], invert=False)

                    right_meshes = working_surface.clip_surface(fault_surfaces[fault][0], invert=True, value=-1*clip_scalar)

                    right_edges = right_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault][0], invert=False, value=-1*(clip_scalar+1))
                    right_edges.compute_normals(cell_normals=True, inplace=True)
                    right_edges_norm = np.average(right_edges['Normals'], axis=0)*flipZ
                        #right_edges

                    print('Right edges avg Norm ' + str(np.around(right_edges_norm, decimals=2)))
                    dirR = (3*clip_scalar*(right_edges_norm + fault_norm)**2)
                    print('Right ' + str(layer) + ' scalar dir ='+ str(np.around(dirR, decimals=2)))
                    right_meshes_ext = right_edges.extrude(dirR).clip_surface(fault_surfaces[fault][0], invert=True)
                        #left_meshes_ext = left_edges.extrude(left_edges['Normals']).clip_surface(fault_surfaces['Main_Fault'][0], invert=True)

                    p.add_mesh(right_meshes_ext, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])
                            #Adding each mesh and extension to the Plotter


                    p.add_mesh(left_meshes_ext, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])
                    working_surface = left_meshes + right_meshes
                    #working_surface += right_meshes 
                    #working_surface += left_meshes
                    #p.add_mesh
                elif fault_norm[0] < 0:
                    print('Fault normal < 0')
                    #fault_norm = faults.compute_normals(cell_normals=True, flip_normals=True)

                    # Average vector of all fault surface normals 

                    #fault_norm = np.average(fault_norm['Normals'], axis=0)
                    left_meshes = working_surface.clip_surface(fault_surfaces[fault][0], invert=True, value=-1*clip_scalar)
                    #left_meshes.compute_normals(cell_normals=False, inplace=True)
                            #finding left edges
                    left_edges = left_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault][0], invert=False, value=-1*(clip_scalar+1))
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
                        left_meshes_ext = left_edges.extrude(directionL).clip_surface(fault_surfaces[fault][0], invert=True)



                        right_meshes = working_surface.clip_surface(fault_surfaces[fault][0], invert=False, value=1*clip_scalar)
                        #right_meshes.compute_normals(cell_normals=False, inplace=True)

                        right_edges = right_meshes.extract_feature_edges().clip_surface(fault_surfaces[fault][0], invert=True, value=1*(clip_scalar+1))

                        #print(right_edgeright_edges = right_meshes.compute_normals(cell_normals=False)s['Normals'])
                        right_edges_norm = np.average(right_edges['Normals'], axis=0)
                            #right_edges

                        print('right edges avg normal' + str(np.around(right_edges_norm, decimals=2)))
                        directionR = (3*clip_scalar*(right_edges_norm - fault_norm)**2)
                        print('Right ' + str(layer) + ' scalar direction ='+ str(np.around(directionR, decimals=2)))
                        right_meshes_ext = right_edges.extrude(directionR).clip_surface(fault_surfaces[fault][0], invert=False)
                            #left_meshes_ext = left_edges.extrude(left_edges['Normals']).clip_surface(fault_surfaces['Main_Fault'][0], invert=True)
                        working_surface = left_meshes + right_meshes
                        p.add_mesh(right_meshes_ext, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])
                                #Adding each mesh and extension to the Plotter


                        p.add_mesh(left_meshes_ext, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])
                        #print('joning meshes')
                    #working_surface += left_meshes
                    #working_surface += right_meshes
            #p.add_mesh(right_meshes, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])          
            #p.add_mesh(left_meshes, color='blue'.to_list()[0])
            p.add_mesh(working_surface, color=modeldf[modeldf['surface']==layer]['color'].to_list()[0])
            for surface in fault_surfaces:
                p.add_mesh(fault_surfaces[surface][0], color=modeldf[modeldf['surface']==surface]['color'].to_list()[0])

    p.show_grid(color='black')
    p.set_background(color='white')
    p.set_scale(1,1,1)

    p.show()
