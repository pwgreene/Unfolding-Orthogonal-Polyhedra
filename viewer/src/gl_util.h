#ifndef GLUTIL_H
#define GLUTIL_H

#include <vector>
#include <vecmath.h>
#include <cstdint>

/* GLProgram is a helper for updating uniform variables.
   Before drawing geometry, update the model matrix and diffuse color.

   You don't have to update the lighting uniforms (they are set at the
   beginning of the frame for you)
*/
class Camera;
struct GLProgram {
    // constructor
    GLProgram(uint32_t program_light, uint32_t program_color, Camera* camera);

    // Update the model matrix. View and projection matrix
    // are read from the camera.
	void updateModelMatrix(Matrix4f M) const;

    // Update material properties.
    // - The one argument version just sets the diffuse color
    // - With 2-3 arguments, also sets specular color
	void updateMaterial(Vector3f diffuseColor, 
        Vector3f ambientColor = Vector3f(-1, -1, -1),
        Vector3f specularColor = Vector3f(0, 0, 0), 
        float shininess = 1.0f,
        float alpha = 1.0f) const;

    // Update lighting. Sets position and color of a single light source
    // in world space.
	void updateLight(Vector3f pos, Vector3f color = Vector3f(1, 1, 1)) const;

    void enableLighting();
    void disableLighting();

private:
    // member variables
    uint32_t active_program;
    uint32_t program_light;
    uint32_t program_color;
    const Camera* camera;
};
#endif
