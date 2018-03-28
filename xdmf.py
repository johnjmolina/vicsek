#!/usr/bin/python
# Python script to generate xdmf file.
# V2.0

# <?xml version="1.0" ?>
# <!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
# <Xdmf Version="2.0">
# <Domain>
#   <Grid GridName="AnyName" GridType="Collection" CollectionType="Temporal">
#
#     <Time TimeType="Hyperslab">
#        <DataItem Dimensions="3" Format="XML" NumberType="Float" Precision="4">
#          t_start t_delta num_t
#        </DataItem>
#     </Time>
#
#     <Grid GridType="Uniform" Name="Particle Frame">
#
#       <!-- Example for unstructured data data -->
#       <Topology NodesPerElement="NUMP" TopologyType="Polyvertex"/>
#       <Geometry GeometryType="XYZ"|"XY">
#         <DataItem Dimensions="NUMP DIM" Format="HDF" NumberType="Float" Precision="4">
#           h5:loc
#         </DataItem>
#       </Geometry>
#
#       <!-- Example for structured data -->
#       <Topology NumberOfElements="Nx Ny Nz" TopologyType="3DSMesh"> 
#       <Geometry GeometryType="X_Y_Z">
#         <DataItem Dimensions="Nx Ny Nz" Format="HDF" NumberType="Float" Precision="4">
#           h5:loc
#         </DataItem>
#         <DataItem Dimensions="Nx Ny Nz" Format="HDF" NumberType="Float" Precision="4">
#           h5:loc
#         </DataItem>
#         <DataItem Dimensions="Nx Ny Nz" Format="HDF" NumberType="Float" Precision="4">
#           h5:loc
#         </DataItem>
#       </Geometry>
#
#       <!-- Example for Scalar Data -->
#       <Attribute AttributeType="Scalar"
#                  Center="Node"|"Edge"|"Face"|"Cell"|"Grid"
#                  Name="ANYNAME">
#         <DataItem Dimensions="NUMP DIM" Format="HDF" NumberType="Float" Precision="4">
#           h5:loc
#         </DataItem>
#       <Attribue>
#
#       <!-- Example for 2d Vector : Paraview requires vector data to be 3D -->
#       <Attribute AttributeType="Vector"  Center="Node" Name="ANYNAME">
#         <DataItem Dimensions="NUMP DIM" ItemType="Function" Function="JOIN($0, $1, 0.0 * $1)">
#           
#           <DataItem Dimension="NUMP" ItemType="HyperSlab" Type="HyperSlab">
#             <DataItem Dimensions="3 DIM" Format="XML"> start_0 start_1 ... stride_0 stride_1... count_0 count_1 />
#             <DataItem Dimensions="NUMP DIM" Format="HDF" Name="Points" NumberType="Float" Precision="4">
#               h5:loc
#             </DataItem>
#           </DataItem>
#
#           <DataItem Dimension="NUMP" ItemType="HyperSlab" Type="HyperSlab">
#             <DataItem Dimensions="3 DIM" Format="XML"> start_0 start_1 ... stride_0 stride_1... count_0 count_1 />
#             <DataItem Dimensions="NUMP DIM" Format="HDF" Name="Points" NumberType="Float" Precision="4">
#               h5:loc
#             </DataItem>
#           </DataItem>
#
#         </DataItem>
#       <Attribute>
#
#       <!-- Example for 3d Vector : -->
#       <Attribute AttributeType="Vector" Center="Node" Name="ANYNAME">
#         <DataItem Dimensions="Nx Ny Nz DIM" ItemType="Function" Function="JOIN($0, $1, $2)">
#           
#           <DataItem Dimension="NUMP" ItemType="HyperSlab" Type="HyperSlab">
#             <DataItem Dimensions="3 DIM" Format="XML"> start_0 start_1 ... stride_0 stride_1... count_0 count_1 />
#             <DataItem Dimensions="NUMP DIM" Format="HDF" Name="Points" NumberType="Float" Precision="4">
#               h5:loc
#             </DataItem>
#           </DataItem>
#
#           <DataItem Dimension="NUMP" ItemType="HyperSlab" Type="HyperSlab">
#             <DataItem Dimensions="3 DIM" Format="XML"> start_0 start_1 ... stride_0 stride_1... count_0 count_1 />
#             <DataItem Dimensions="NUMP DIM" Format="HDF" Name="Points" NumberType="Float" Precision="4">
#               h5:loc
#             </DataItem>
#           </DataItem>
#
#         </DataItem>
#       <Attribute>
#
#     </Grid>
#
#   </Grid>
# </Domain>
# </Xdmf>
import xml.etree.ElementTree as et
from xml.dom import minidom

#valid types: float, int, uint, char, uchar
real4 = {'type':'Float', 'prec':'4'}
real8 = {'type':'Float', 'prec':'8'}
int4  = {'type':'Int',   'prec':'4'}
int8  = {'type':'Int',   'prec':'8'}

#tensor types: Scalar(1), Vector(3), Tensor(9), Tensor(6)
scalar = {'type':'Scalar', 'rank':1}
vector = {'type':'Vector', 'rank':3, 'join':'JOIN($0, $1, $2)', 'join2d':'JOIN($0, $1, (0.0 * $1))'}
tensor = {'type':'Tensor', 'rank':9, 'join':'JOIN($0, $1, $2, $3, $4, $5, $6, $7, $8)'}
stensor= {'type':'Tensor6','rank':6, 'join':'JOIN($0, $1, $2, $3, $4, $5)'}



def isIter(l):
    return isinstance(l, (list, tuple))

def printList(l):
    if isIter(l):
        s=str(l[0])
        for i in range(1, len(l)):
            s+=" "+str(l[i])
        return s
    else:
        return str(l)
    
def addElem(elem, name, params={}, text=None):
    elem.append(et.Element(name))
    for key,val in params.iteritems():
        elem[-1].set(key,val)
    if text is not None:
        elem[-1].text = text
    return [elem[-1]]

def addSubElem(elem, name, params={}, text=None):
    elem.append(et.SubElement(elem[-1], name))
    for key,val in params.iteritems():
        elem[-1].set(key,val)
    if text is not None:
        elem[-1].text = text
    return [elem[-1]]

def addDataItem(elem, ndims, dtype=real4, format='HDF', text=None):
    if text is not None:
        addSubElem(elem, 'DataItem',
                   params={'Dimensions':ndims, 'Format':format,
                           'NumberType':dtype['type'], 'Precision':dtype['prec']},
                   text=text)
    else:
        addSubElem(elem, 'DataItem',
                   params={'Dimensions':ndims, 'Format':format,
                           'NumberType':dtype['type'], 'Precision':dtype['prec']},
                   text=text)
def addDataSlab(elem, grid_dims, slab_dims, data_dims, slab_text, text, dtype=real4):
    addSubElem(elem, 'DataItem', params={'Dimensions':grid_dims, 'ItemType':'Hyperslab', 'Type':'HyperSlab'})
    addSubElem(elem, 'DataItem', params={'Dimensions':slab_dims, 'Format':'XML'}, \
               text=slab_text)
    elem.pop()
    addDataItem(elem, data_dims, dtype=dtype, text=text)
    elem.pop()
    elem.pop()
    
def addFunction(elem, ndims, func):
    addSubElem(elem, 'DataItem',\
               params={'Dimensions':ndims, 'ItemType':'Function', 'Function':func})
    
        
def free(elem, fname):
    """ Finalize and xdmf xml tree 
    Args: 
        elem (xml)
    """
    
    def prettyprintxml(elem):
        rough_string = et.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        reparsed = reparsed.toprettyxml(indent="  ")
        return reparsed

    rough_string = "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" + et.tostring(elem[0], 'utf-8')
    reparsed = minidom.parseString(rough_string)
    fout = open(fname, 'w')
    fout.write(reparsed.toprettyxml(indent="  "))
    fout.close()

def init(time, name):
    """ Initialize an xdmf xml tree
    Args: 
        time ([string, string, string]): array with t_start, t_delta, and num_frames for time series data
        name (string): grid data name
    Returns:
        xdmf: bare xml tree
    """
    xdmf = []
    addElem(xdmf, 'Xdmf', params={'Version':'2.0'})
    addSubElem(xdmf, 'Domain')
    grid = addSubElem(xdmf, 'Grid', params={'CollectionType':'Temporal', 'GridType':'Collection', 'Name':name})
    addSubElem(grid, 'Time', params={'TimeType':'HyperSlab'})
    addDataItem(grid, '3', dtype=real4, format='XML', text=printList(time))

    return xdmf

def __newGridND(xdmf, conf, nloc, name, dtype=real4):
    grid  = addSubElem(xdmf, 'Grid', params={'GridType':'Uniform', 'Name':name})    
    addSubElem(grid, 'Topology', params={'NumberOfElements':conf['grid_dims'], 'TopologyType':conf['topology']})
    grid.pop()
    addSubElem(grid, 'Geometry', params={'GeometryType':conf['geometry']})
    for loc in nloc:
        addDataItem(grid, conf['data_dims'], dtype=dtype, text=loc)
        grid.pop()
    return grid

def newGrid3D(xdmf, ns, nloc, name, dtype=real4):
    assert(isIter(ns) and isIter(nloc) and len(ns) == 3 and len(nloc) == 3)
    dims = printList(ns)    
    conf = {'topology':'3DSMesh', 'geometry':'X_Y_Z',
            'grid_dims':dims, 'data_dims':dims,
            'slab_dims':'3 4',
            'slab_x':printList([0, 0, 0, 0, 1, 1, 1, 1, ns[0], ns[1], ns[2], 1]),
            'slab_y':printList([0, 0, 0, 1, 1, 1, 1, 1, ns[0], ns[1], ns[2], 1]),
            'slab_z':printList([0, 0, 0, 2, 1, 1, 1, 1, ns[0], ns[1], ns[2], 1])}
    grid = __newGridND(xdmf, conf, nloc, name, dtype=dtype)
    return conf

def newGrid2Das3D(xdmf, ns, nloc, name, dtype=real4):
    assert(isIter(ns) and isIter(nloc) and len(ns) == 2 and len(nloc) == 2)
    dims = printList(ns)
    conf = {'topology':'3DSMesh', 'geometry':'X_Y_Z',
            'grid_dims':dims+' 1', 'data_dims':dims,
            'slab_dims':'3 3',
            'slab_x':printList([0, 0, 0, 1, 1, 1, ns[0], ns[1], 1]),
            'slab_y':printList([0, 0, 1, 1, 1, 1, ns[0], ns[1], 1])}    
    grid = __newGridND(xdmf, conf, nloc, name, dtype=dtype)
    addFunction(grid, conf['data_dims'], '(0.0 * $0)')
    addDataItem(grid, conf['data_dims'], dtype=dtype, text=nloc[-1])
    return conf

def newGrid2D(xdmf, ns, nloc, name, dtype=real4):
    assert(isIter(ns) and isIter(nloc) and len(ns) == 2 and len(nloc) == 2)
    dims = printList(ns)
    conf = {'topology':'2DSMesh', 'geometry':'X_Y',
            'grid_dims':dims+' 1', 'data_dims':dims,
            'slab_dims':'3 3',
            'slab_x':printList([0, 0, 0, 1, 1, 1, ns[0], ns[1], 1]),
            'slab_y':printList([0, 0, 1, 1, 1, 1, ns[0], ns[1], 1])}    
    grid = __newGridND(xdmf, conf, nloc, name, dtype=dtype)
    return conf

def __newPolyND(xdmf, conf, nloc, name, dtype=real4):
    grid = addSubElem(xdmf, 'Grid', params={'GridType':'Uniform', 'Name':name})
    addSubElem(grid, 'Topology', \
               params={'NodesPerElement':conf['poly_dims'],
                       'TopologyType':conf['topology']})
    grid.pop()
    addSubElem(grid, 'Geometry', params={'GeometryType':conf['geometry']})
    if isinstance(nloc, basestring):
        addDataItem(grid, conf['poly_dims'], dtype=dtype, text=nloc)
    elif isIter(nloc):
        addFunction(grid, conf['poly_dims'], conf['poly_join'])
        for loc in nloc:
            addDataItem(grid, conf['data_dims'], dtype=dtype, text=loc)
            grid.pop()
    return grid

def newPoly3D(xdmf, ns, nloc, name, dtype=real4):
    assert(not isIter(ns) and (not isIter(nloc) or (isIter(nloc) and len(nloc) == 3)))
    dims = printList(ns)
    conf = {'topology':'Polyvertex', 'geometry':'XYZ',
            'poly_dims':dims+' 3', 'grid_dims':dims, 'data_dims':dims, 
            'poly_join':'JOIN($0, $1, $2)',
            'slab_dims':'3 2',
            'slab_x':printList([0, 0, 1, 1, ns, 1]),
            'slab_y':printList([0, 1, 1, 1, ns, 1]),
            'slab_z':printList([0, 2, 1, 1, ns, 1])}
    grid = __newPolyND(xdmf, conf, nloc, name, dtype=dtype)
    return conf

def newPoly2D(xdmf, ns, nloc, name, dtype=real4):
    assert(not isIter(ns) and (not isIter(nloc) or (isIter(nloc) and len(nloc) == 2)))
    dims = printList(ns)
    conf = {'topology':'Polyvertex', 'geometry':'XY',
            'poly_dims':dims+' 2', 'grid_dims':dims, 'data_dims':dims, 
            'poly_join':'JOIN($0, $1)',
            'slab_dims':'3 2',
            'slab_x':printList([0, 0, 1, 1, ns, 1]),
            'slab_y':printList([0, 1, 1, 1, ns, 1])}
    grid = __newPolyND(xdmf, conf, nloc, name, dtype=dtype)
    return conf

def freeGrid(xdmf):
    xdmf.pop()
def freePoly(xdmf):
    xdmf.pop()
            
def addTensorData(grid, conf, nloc, ttype, name, dtype=real4):
    tdims = conf['grid_dims'] + ' ' + str(ttype['rank'])
    addSubElem(grid, 'Attribute', params={'AttributeType':ttype['type'], 'Center':'Node', 'Name':name})
    if isinstance(nloc, basestring):
        addDataItem(grid, tdims, dtype=dtype, text=nloc)
    elif isIter(nloc) and (ttype != scalar) and (len(nloc) == ttype['rank']):
        addFunction(grid, tdims, ttype['join'])
        for loc in nloc:
            addDataItem(grid, conf['data_dims'], dtype=dtype, text=loc)
            grid.pop()
    else:
        print("Invalid Tensor Specification")
        exit(1)
    grid.pop()
    grid.pop()

def addScalar(grid, conf, nloc, name, dtype=real4):
    addTensorData(grid, conf, nloc, scalar, name, dtype)

def addVector(grid, conf, nloc, name, dtype=real4):
    addTensorData(grid, conf, nloc, vector, name, dtype)

def addVector2D(grid, conf, nloc, name, dtype=real4):
    v3dims = conf['grid_dims'] + ' 3'
    v2dims = conf['data_dims'] + ' 2'
    addSubElem(grid, 'Attribute', \
               params={'AttributeType':vector['type'], 'Center':'Node', 'Name':name})
    addFunction(grid, v3dims, vector['join2d'])    
    if isinstance(nloc, basestring):
        for slab in ['slab_x', 'slab_y']:
            addDataSlab(grid, conf['data_dims'], conf['slab_dims'], v2dims, conf[slab], nloc, dtype=dtype)
    elif isIter(nloc) and  len(nloc) == 2:
        for loc in nloc:
            addDataItem(grid, conf['data_dims'], dtype=dtype, text=loc)
            grid.pop()
    else:
        print("Invalid Vector2D Specification")
        exit(1)
    grid.pop()
    grid.pop()

def addSymTensor(grid, conf, nloc, name, dtype=real4):
    addTensorData(grid, conf, nloc, stensor, name, dtype)

def addTensor(grid, conf, nloc, name, dtype=real4):
    addTensorData(grid, conf, nloc, tensor, name, dtype)

def addMatrix(grid, conf, mrank, loc, name, dtype=real4):
    mdims = conf['data_dims'] + ' ' + printList(mrank)
    addSubElem(grid, 'Attribute', params={'AttributeType':'Matrix', 'Center':'Node', 'Name':name})
    addDataItem(grid, mdims, dtype=dtype, text=loc)
    grid.pop()
    grid.pop()
