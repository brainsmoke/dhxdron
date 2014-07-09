
from math import *
from linear import *
import svg

#
# settings
#

radius = 160.
thickness = .9


#
# inkscape sizes
#

dpi = 90.
native_scale = dpi/25.4

radius *= native_scale
thickness *= native_scale

#

jagged_longedge = tuple( (x/10., y) for x,y in (
    (0,  0),
    (2,  0),
    (2,  1.5),
    (3,  1.5),
    (3, -1.5),
    (3.15, -1.5),
    (3.15, -1),
    (3.85, -1),
    (3.85, -1.5),
    (4, -1.5),
    (4,  0),
    (6,  0),
    (6,  1.5),
    (7,  1.5),
    (7, -1.5),
    (7.15, -1.5),
    (7.15, -1),
    (7.85, -1),
    (7.85, -1.5),
    (8, -1.5),
    (8,  0),
    #(10, 0),
) )


slot_long = tuple( (x/10., y) for x,y in (
    (4.8,  -3),
    (4.8,  -4.2),
    (5.2,  -4.2),
    (5.2,  -3),
) )

jagged_shortedge = tuple( (x/7., y) for x,y in (
    (0,  0),
    (1,  0),
    (1,  1.5),
    (2,  1.5),
    (2, -1.5),
    (2.15, -1.5),
    (2.15, -1),
    (2.85, -1),
    (2.85, -1.5),
    (3, -1.5),
    (3,  0),
    (4,  0),
    (4,  1.5),
    (5,  1.5),
    (5, -1.5),
    (5.15, -1.5),
    (5.15, -1),
    (5.85, -1),
    (5.85, -1.5),
    (6, -1.5),
    (6,  0),
    #(7, 0),
) )

slot_short = tuple( (x/7., y) for x,y in (
    (3.3,  -3),
    (3.3,  -4.2),
    (3.7,  -4.2),
    (3.7,  -3),
) )

def normalize2( (x, y) ):
    d = sqrt(x*x+y*y)
    return (x/d, y/d)

def vector_sub2( (x1,y1), (x2,y2) ):
    return (x1-x2, y1-y2)

def vector_add2( (x1,y1), (x2,y2) ):
    return (x1+x2, y1+y2)

def interpolate2( (x1,y1), (x2,y2), frac ):
    return ( x2*frac+x1*(1-frac), y2*frac+y1*(1-frac) )

def jagedge(a, b, jag, w):
    edges = []

    x, y = a
    dx, dy = normalize2(vector_sub2(b, a))
    for l, v in jag:
        edges.append( vector_add2(interpolate2(a, b, l), (-dy*w*v, dx*w*v) ) )

    return edges
        

def scale_to_plane_on_normal(point, normal):
    return scalar_mul(scalar_product(normal,normal) / scalar_product(point,normal), point)

def dhxdron_deltoid_coords():

    phi = ( sqrt(5.) + 1. ) / 2.
    phi2 = phi**2
    phi3 = phi**3

    #rhombicosidodecahedron

    normal = normalize( ( 1,  1, phi3) )

    square1 = (
        ( 1,  1, phi3),
        (-1,  1, phi3),
        ( 1, -1, phi3),
        (-1, -1, phi3),
    )

    triangle = (
        ( 1,  1, phi3),
        (-1,  1, phi3),
        (0, phi2, (2+phi)),
    )

    pentagon = (
        ( 1,  1, phi3),
        (phi2, phi, 2*phi),
        ((2+phi), 0, phi2),
        (phi2, -phi, 2*phi),
        ( 1, -1, phi3),
    )

    square2 = (
        ( 1,  1, phi3),
        (0, phi2, (2+phi)),
        (phi, 2*phi, phi2),
        (phi2, phi, 2*phi),
    )

    deltoid = (
        scale_to_plane_on_normal(vector_sum( *triangle ), normal),
        scale_to_plane_on_normal(vector_sum( *square1 ), normal),
        scale_to_plane_on_normal(vector_sum( *pentagon ), normal),
        scale_to_plane_on_normal(vector_sum( *square2 ), normal),
    )

    return deltoid

def dhxdron_deltoid_2d_coords():
    deltoid = dhxdron_deltoid_coords()

    short_axis   = d(deltoid[1], deltoid[3])
    long_axis    = d(deltoid[0], deltoid[2])
    cross        = interpolate(deltoid[1], deltoid[3], .5)
    sym_ax_long  = d(cross, deltoid[2])
    sym_ax_short = d(cross, deltoid[0])

    return ( (             0.,  long_axis/2.              ),
             (  short_axis/2.,  long_axis/2.-sym_ax_short ),
             (             0., -long_axis/2.              ),
             ( -short_axis/2.,  long_axis/2.-sym_ax_short ) )

def deltoid_get_angle():
    deltoid = dhxdron_deltoid_coords()

    short_axis   = d(deltoid[1], deltoid[3])
    long_axis    = d(deltoid[0], deltoid[2])
    cross        = interpolate(deltoid[1], deltoid[3], .5)
    sym_ax_long  = d(cross, deltoid[2])
    sym_ax_short = d(cross, deltoid[0])

    x, y, z = deltoid[3]
    long_diag = d(deltoid[1], deltoid[2])
    short_diag = d(deltoid[0], deltoid[1])

    shortcut_long = y
    ortholine_long = sym_ax_long*short_axis/long_diag
    return acos( shortcut_long/ortholine_long )*2

    # SAME! :-)
    #shortcut_short = x
    #ortholine_short = sym_ax_short*short_axis/short_diag
    #print acos(shortcut_short/ortholine_short)*2

def get_angle():
    """ alternative calculation for sanity check """
    phi = ( sqrt(5.) + 1. ) / 2.

    r = d( (0,0,0), (1,1,phi**3) )
    return acos(  ( r**2 + r**2 - 2**2 ) / ( r*r*2 ) )


def shape(radius, thickness):
    edges = []
    a, b, c, d = ( (x*radius, y*radius) for x,y in dhxdron_deltoid_2d_coords() )
    
    angle = deltoid_get_angle()
    indent = 1./tan(angle)*thickness

    edges = (
        jagedge(a, b, jagged_shortedge, indent),
        jagedge(b, c, jagged_longedge, indent),
        jagedge(c, d, jagged_longedge, indent),
        jagedge(d, a, jagged_shortedge, indent),
    )
    return [ c for e in edges for c in e ]


def slots(radius):
    edges = []
    a, b, c, d = ( (x*radius, y*radius) for x,y in dhxdron_deltoid_2d_coords() )
    
    return (
        jagedge(a, b, slot_short, native_scale),
        jagedge(b, c, slot_long, native_scale),
        jagedge(c, d, slot_long, native_scale),
        jagedge(d, a, slot_short, native_scale),
    )

style = 'stroke:none;fill:#0000ff;opacity:.3'

print svg.header(radius*2,radius*2)

print '<g transform="translate('+str(radius)+' '+str(radius)+')">'
print svg.circle(radius, style)
print svg.path( svg.polygon_path(shape(radius=radius, thickness=thickness)), style)
for slot in slots(radius=radius):
    print svg.path( svg.polygon_path(slot), style)

print '</g>'
print svg.footer()

