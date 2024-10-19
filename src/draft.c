#include <surface/surface.h>
#include <stdlib.h>

// let's make a cube out of 4 squares


// note that a square is a function 

        // it would be nice to have parameterizable functions
        // like to write a function that takes 5 parameters,
        // and set 3 of them, and treat it as a surface

        // it would be nice to do it small,
        // ie to either have lambdas or local functions,
        // or to pass the 3 parameters to a further function to
        // get a pointer to a function.

            // it looks like that the only reasonable way to do that
            // is to make a global-scoped function for each use.
            // this function could reference data from a global or passed object; surfaces could take userdata to produce that

typedef oriented_surface square3;
pt3 square3_func(surface * surface, pt2 uv)
{
    square3 * square = (square3*)surface;
    pt3 pt = square->position;
    for (int dim = 0; dim < 2; ++ dim)
        add3into(&pt, pt, mul31(square->axes[dim], uv.idx[dim] * 2 - 1));
    return pt;
}
square3 square3_new(pt3 position, pt3 axis1, pt3 axis2)
{
    return (square3){
        .surface = {.func = square3_func, .neighbors = {0}},
        .position = position,
        .axes = {axis1, axis2}
    };
}

typedef surface_group cube;
void cube_free(surface_group* cube)
{
    free(*cube->surfaces);
    free(cube->surfaces);
}
cube cube_new(pt3 position, quat orientation_and_size)
{
    oriented_surface * surfaces = malloc(sizeof(oriented_surface)*6);
    surface ** surface_list = malloc(sizeof(surface*)*6);
    pt3 axes[3];
    transformationaxesinto(&axes, orientation_and_size);
    // iterate over 6 faces with idcs
    // each face has a sign pair, and a direction.
    const int FACE_AXIS_STRIDE = 1;
    const int FACE_SIGN_STRIDE = 3;
    for (int sign = 1; sign > -2; sign -= 2) {
        for (int axis = 0; axis < 3; ++ axis) {
            int idx = axis*FACE_AXIS_STRIDE+(sign+1)*FACE_SIGN_STRIDE/2;
            surface_list[idx] = &surfaces[idx].surface;
            surfaces[idx] = square3_new(
                add3(position, mul31(axes[axis], sign)), 
                axes[(axis+1)%3],
                axes[(axis+2)%3]
            );
            for (int neighbor_axis = 0; neighbor_axis < 2; ++ neighbor_axis) {
                for (int neighbor_dir = 0; neighbor_dir < 2; ++ neighbor_dir) {
                    int neighbor_idx = (
                        ((axis+neighbor_axis)%3) * FACE_AXIS_STRIDE +
                        neighbor_dir * FACE_SIGN_STRIDE
                    );
                    *surface_neighbor(&surfaces[idx].surface, neighbor_axis, neighbor_dir) = &surfaces[neighbor_idx].surface;
                }
            }
        }
    }
    return (cube){
        .surfaces = surface_list,
        .n_surfaces = 6,
        .free = cube_free,
    };
}

int main()
{
    cube a_cube = cube_new((pt3){0,0,0},(quat){0,1,0,0});
}
