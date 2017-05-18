#include "gl.h"
#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstdio>
#include <cstdint>

#include "vertexrecorder.h"
#include "recorder.h"
#include "starter3_util.h"
#include "camera.h"
#include "gl_util.h"

using namespace std;

namespace
{
    
    // Declarations of functions whose implementations occur later.

    //void initRendering();

    // Some constants
    const Vector3f LIGHT_POS(0.0f, 8.0f, 16.0f);
    const Vector3f LIGHT_COLOR(500, 500, 500);
    const Vector3f MESH_COLOR(0.5f, 0.5f, 0.5f);
    const Vector3f MESH_AMBIENT = .3 * MESH_COLOR;
    Vector3f PATH_COLOR(0.0f, 0.0f, 0.0f);
    Vector3f PATH_AMBIENT = .9 * PATH_COLOR;
    const float ambient_amount = .15;
    
    bool MESH_ON = false;

    
    // Globals here.

    // for rendering a mesh
    vector<Vector3f> vecv;
    vector<Vector3f> vecn;
    vector<vector<unsigned>> vecf;
    
    //for rendering path(s)
    vector<vector<Vector3f>> paths;
    
    float PATH_WIDTH = .01;
    float SCALE = 10;

    Camera camera;
    bool gMousePressed = false;
    GLuint program_color;
    GLuint program_light;

    // Function implementations
    static void keyCallback(GLFWwindow* window, int key,
        int scancode, int action, int mods)
    {
        if (action == GLFW_RELEASE) { // only handle PRESS and REPEAT
            return;
        }

        // Special keys (arrows, CTRL, ...) are documented
        // here: http://www.glfw.org/docs/latest/group__keys.html
        switch (key) {
        case GLFW_KEY_ESCAPE: // Escape key
            exit(0);
            break;
        case ' ':
        {
            Matrix4f eye = Matrix4f::identity();
            camera.SetRotation(eye);
            camera.SetCenter(Vector3f(0, 0, 0));
            break;
        }
        case 77:
            MESH_ON = !MESH_ON; //toggle mesh appearance
            break;
        default:
            cout << "Unhandled key press " << key << "." << endl;
        }
    }

    static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
    {
        double xd, yd;
        glfwGetCursorPos(window, &xd, &yd);
        int x = (int)xd;
        int y = (int)yd;

        int lstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        int rstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
        int mstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE);
        if (lstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::LEFT, x, y);
        }
        else if (rstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::RIGHT, x, y);
        }
        else if (mstate == GLFW_PRESS) {
            gMousePressed = true;
            camera.MouseClick(Camera::MIDDLE, x, y);
        }
        else {
            gMousePressed = true;
            camera.MouseRelease(x, y);
            gMousePressed = false;
        }
    }

    static void motionCallback(GLFWwindow* window, double x, double y)
    {
        if (!gMousePressed) {
            return;
        }
        camera.MouseDrag((int)x, (int)y);
    }

    void setViewport(GLFWwindow* window)
    {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);

        camera.SetDimensions(w, h);
        camera.SetViewport(0, 0, w, h);
        camera.ApplyViewport();
    }

    
    void loadMesh()
    {
        //for .off files
        const int MAX_BUFFERSIZE = 4096;
        char buffer[MAX_BUFFERSIZE];
        cin.getline(buffer, MAX_BUFFERSIZE);
        cin.getline(buffer, MAX_BUFFERSIZE);
        stringstream ssn(buffer);
        
        vector<vector<int>> adjacentFaces; //keep track of adjacent faces to each vertex
        int nv, nt, ne;
        vector<Vector3f> faceNormals;
        ssn >> nv >> nt >> ne;
        printf("%d, %d, %d\n", nv, nt, ne);
        adjacentFaces.resize(nv);
        for (int v_i = 0; v_i < nv; v_i++) {
            Vector3f v;
            Vector3f n;
            cin.getline(buffer, MAX_BUFFERSIZE);
            stringstream ss(buffer);
            ss >> v[0] >> v[1] >> v[2];
            vecv.push_back(v);
        }
        printf("loaded %d vertices\n", nv);
        for (int f_i = 0; f_i < nt; f_i++) {
            vector<unsigned> f;
            unsigned x, a, b, c;
            cin.getline(buffer, MAX_BUFFERSIZE);
            stringstream ss(buffer);
            ss >> x >> a >> b >> c; //just assuming triangle mesh for simplicity
            f.push_back(a); f.push_back(a); f.push_back(a);
            f.push_back(b); f.push_back(b); f.push_back(b);
            f.push_back(c); f.push_back(c); f.push_back(c);
            if (a > nv || b > nv || c > nv) {
                printf("bad, %d %d %d\n", a, b, c);
            }
            vecf.push_back(f);
            //compute face normals
            Vector3f ij = vecv[f[3]] - vecv[f[0]];
            Vector3f jk = vecv[f[6]] - vecv[f[3]];
//            X[f[0]] = ij.normalized(); X[f[3]] = -ij.normalized(); X[f[6]] = -jk.normalized(); //arbitrarily choose these--for now
            Vector3f n = Vector3f::cross(ij, jk);
            faceNormals.push_back(n.normalized());
            for (int i = 0; i < 9; i+=3) {
                adjacentFaces[f[i]].push_back(f_i);
            }
        }
        printf("loaded %d faces\n", nt);
        for (int v_i = 0; v_i < nv; v_i++) {
            Vector3f v_n = Vector3f(0, 0, 0);
            
            int nfaces = adjacentFaces[v_i].size();
            if (nfaces == 0) {
                printf("%d\n", v_i);
            }
            for (int i = 0; i < nfaces; i++) {
                v_n += faceNormals[adjacentFaces[v_i][i]];
            }
            v_n /= nfaces;
            v_n.normalize();
            
            vecn.push_back(v_n);
        }
    }
    
    void loadPath() {
        const int MAX_BUFFERSIZE = 4096;
        char buffer[MAX_BUFFERSIZE];
        cin.getline(buffer, MAX_BUFFERSIZE);
        stringstream ss(buffer);
        
        int nPaths, pathLen;
        float x, y, z;
        ss >> nPaths;
        
        //load each path one at a time from stdin
        for (int i = 0; i < nPaths; i++) {
            cin.getline(buffer, MAX_BUFFERSIZE);
            stringstream ss(buffer);
            ss >> pathLen;
            vector<Vector3f> path;
            
            for (int v = 0; v < pathLen; v++) {
                cin.getline(buffer, MAX_BUFFERSIZE);
                stringstream ssv(buffer);
                ssv >> x >> y >> z;
                Vector3f vertex(x, y, z);
                path.push_back(vertex);
            }
            paths.push_back(path);
        }
        
        printf("loaded %d paths\n", nPaths);
        
//        prints out all paths, for testing
        //for (int i = 0; i < nPaths; i++) {
         //   printf("%d\n", paths[i].size());
           // for (int j = 0; j < paths[i].size(); j++) {
             //   paths[i][j].print();
           // }
            //printf("\n");
        }
    }
    
    bool approxEqual(float a, float b) {
        float eps = 1e-4;
        return (a < b + eps && a > b - eps);
    }
    
    void drawPath() {
        GLProgram gl(program_light, program_color, &camera);
        gl.updateMaterial(PATH_COLOR, PATH_AMBIENT);
        Vector3f a(0, 1, 0); //default orientation of cylinder
        for (int i = 0; i < paths.size(); i++) {
            for (int j = 0; j < paths[i].size() - 1; j++) {
                
                Vector3f direction = paths[i][j+1] - paths[i][j]; //direction of path segment
                float length = direction.abs();
                direction.normalize();
                Vector3f location = paths[i][j];
                Vector3f axis = Vector3f::cross(a, direction);
                float dot = Vector3f::dot(a, direction);
                
                //rotates cylinder to correct direction
                Matrix4f rotation;
                if (approxEqual(dot, 1)) {//already aligned with +y-axis
                    rotation = Matrix4f::identity();
                } else if (approxEqual(dot, -1)) { //need to flip
                    rotation = Matrix4f::rotateX(M_PI);
                } else {
                    float rotationAngle = acos(dot);
                    rotation = Matrix4f::rotation(axis, rotationAngle);
                }
                
                //translates cylinder to correct position
                Matrix4f translation = Matrix4f::translation(location);
                //if scaling:
                //translation = scale*translation
                
                gl.updateModelMatrix(translation*rotation);
                drawCylinder(3, PATH_WIDTH, length);
            }
        }
    }
    
    void drawMesh() {
        // draw obj mesh here
        // read vertices and face indices from vecv, vecn, vecf
        GLProgram gl(program_light, program_color, &camera);
        gl.updateModelMatrix(Matrix4f::identity());
        gl.updateLight(LIGHT_POS, LIGHT_COLOR.xyz()); // once per frame
        gl.updateMaterial(MESH_COLOR, MESH_AMBIENT);
        GeometryRecorder rec(vecf.size() * 3);
        for (vector<unsigned> f : vecf) {
            rec.record(Vector3f(vecv[f[0]][0], vecv[f[0]][1], vecv[f[0]][2]),
                       Vector3f(vecn[f[2]][0], vecn[f[2]][1], vecn[f[2]][2]));
            
            rec.record(Vector3f(vecv[f[3]][0], vecv[f[3]][1], vecv[f[3]][2]),
                       Vector3f(vecn[f[5]][0], vecn[f[5]][1], vecn[f[5]][2]));
            
            rec.record(Vector3f(vecv[f[6]][0], vecv[f[6]][1], vecv[f[6]][2]),
                       Vector3f(vecn[f[8]][0], vecn[f[8]][1], vecn[f[8]][2]));
        }
        rec.draw();
        
    }
        
        

    //-------------------------------------------------------------------

    void initRendering()
    {
        // Clear to black
        glClearColor(0, 0, 0, 1);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    // Main routine.
    // Set up OpenGL, define the callbacks and start the main loop
    int main(int argc, char** argv)
    {

        GLFWwindow* window = createOpenGLWindow(1024, 1024, "Assignment 3");

        // setup the event handlers
        glfwSetKeyCallback(window, keyCallback);
        glfwSetMouseButtonCallback(window, mouseCallback);
        glfwSetCursorPosCallback(window, motionCallback);

        initRendering();

        // The program object controls the programmable parts
        // of OpenGL. All OpenGL programs define a vertex shader
        // and a fragment shader.
        program_color = compileProgram(c_vertexshader, c_fragmentshader_color);
        if (!program_color) {
            printf("Cannot compile program\n");
            return -1;
        }
        program_light = compileProgram(c_vertexshader, c_fragmentshader_light);
        if (!program_light) {
            printf("Cannot compile program\n");
            return -1;
        }

        camera.SetDimensions(600, 600);
        camera.SetPerspective(50);
        camera.SetDistance(20);

        loadMesh();
        loadPath();
        
        // Main Loop
        while (!glfwWindowShouldClose(window)) {
            // Clear the rendering window
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            setViewport(window);
            
            if (MESH_ON) {
                PATH_COLOR = Vector3f(0.9f, 0.9f, 0.9f);
                PATH_AMBIENT = .9 * PATH_COLOR;
                drawMesh();
            } else {
                PATH_COLOR = Vector3f(0.0f, 0.0f, 0.0f);
                PATH_AMBIENT = .9 * PATH_COLOR;
            }
            
            drawPath();

            // Make back buffer visible
            glfwSwapBuffers(window);

            // Check if any input happened during the last frame
            glfwPollEvents();
        }

        // All OpenGL resource that are created with
        // glGen* or glCreate* must be freed.
        glDeleteProgram(program_color);
        glDeleteProgram(program_light);


        return 0;	// This line is never reached.
    }
