#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <learnopengl/shader_m.h>
#include <learnopengl/camera.h>
#include <learnopengl/model.h>

#include <iostream>

#include "structures.h"
#include "tests.h"
#include "functions.h"
#include "obj_parser.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

bool WireframeMode = false;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // tell stb_image.h to flip loaded texture's on the y-axis (before loading model).
    // stbi_set_flip_vertically_on_load(true);

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile shaders
    // -------------------------
    Shader ourShader("C:\\Users\\Arcasha\\Desktop\\tmp\\cut\\Clip\\shaders\\model.vs", "C:\\Users\\Arcasha\\Desktop\\tmp\\cut\\Clip\\shaders\\model.fs", WireframeMode);

    // INITIAL TRIANGULATED MODEL
    mdl::Model ourModel("C:\\Users\\Arcasha\\Desktop\\tmp\\cut\\Clip\\models\\cube.obj");
    //mdl::Model ourModel2("C:\\Users\\Arcasha\\Desktop\\tmp\\cut\\Clip\\models\\new.obj");

    // ...........................................................................................................s....................................................

    Mesh m = Convert(ourModel.meshes[0]);

    //Mesh m2 = Convert(ourModel2.meshes[0]);


    // Vector r = {0.5,0.5,0.5};
    // double s = 1.5;
    // Scalation(m, s,s, s);
    // Translation(m2,r);



    Test(m, precise);

    Flat plane2;


    // Dodecahedron( m, precise);

    // Icosahedron( m, precise);

    // Cuboctahedron( m, precise );

    // plane2.p = {{0}, {0}, {1 - precise*9}};//not change
    // plane2.n = Vector{{0}, {0}, {1}}.normalize();
    // Intersect(m, plane2, precise);

    // plane2.p = {{0}, {1}, {1 - precise*17}};//not change
    // plane2.n = Vector{{1}, {1}, {1}}.normalize();
    // Intersect(m, plane2, precise);

    // plane2.p = {{0}, {1 - precise * 9}, {1}};//not change
    // plane2.n = Vector{{0}, {1}, {0}}.normalize();
    // Intersect(m, plane2, precise);

    // double scalar = 400000;//при 1000 уже изменяется кол-во вершин и граней. при 50 уже некорректно. при 10 уже Empty intersect
    // plane2.p = {{0}, {0}, {1 - precise*scalar}};// change
    // plane2.n = Vector{{0}, {0}, {-1}}.normalize();
    // Intersect(m, plane2, precise);

    // plane2.p = {{0}, {1}, {1 - precise*scalar}};// change
    // plane2.n = Vector{{-1}, {-1}, {-1}}.normalize();
    // Intersect(m, plane2, precise);

    // plane2.p = {{0}, {1 - precise * scalar}, {1}};// change
    // plane2.n = Vector{{0}, {-1}, {0}}.normalize();
    // Intersect(m, plane2, precise);




    // Rhombicuboctahedron( m, precise );

    // Rhombicuboctahedron2( m, precise );

    // Rhombicuboctahedron3( m, precise );

    // Pyramid( m, precise );

    // Octahedron(m, precise);

    // Tetrahedron(m, precise);




    plane2.p = {{0}, {0}, {0}};
    plane2.n = Vector{{1}, {1}, {1}}.normalize();
    Intersect(m, plane2, precise);


    // plane2.p = {{0}, {0}, {0}}; 
    // plane2.n = Vector{{1}, {0}, {0}}.normalize();
    // Intersect(m, plane2, precise);

    // plane2.p = {{0}, {0.7}, {0}}; 
    // plane2.n = Vector{{0}, {1}, {0}}.normalize();
    // Intersect(m, plane2, precise);




    Write(m, precise);//запись в new.obj полученной модели


    std::vector<mdl::Vertex> ver; //convert Mesh back to mdl::Mesh
    std::vector<Vertex> m_ver;

    for(int i = 0; i < m.Faces.size(); ++i){//пересчитать нормали
        Vector n = Normal( m.Faces[i], m, precise);
        for(int j = 0; j < m.Faces[i].Indices.size(); ++j){
            Vertex v = m.Vertices[m.Faces[i].Indices[j]];
            mdl::Vertex vert;//mdl
            vert.Position = {v.x, v.y, v.z};
            vert.Normal = {n.x, n.y, n.z};
            m_ver.push_back(v);
            ver.push_back(vert);
            m.Faces[i].Indices[j] = m_ver.size() - 1;
        }
    }
    std::vector<unsigned int> indices;
    for(int i =0; i < m.Faces.size(); ++i)
        for( int j = 0; j < m.Faces[i].Indices.size(); ++j)
            indices.push_back(m.Faces[i].Indices[j]);

    bool IsTr = Check(m);
    std::cout<<"v: "<<m.Vertices.size() <<"   f: "<<m.Faces.size()<<std::endl;
    std::cout<<"is correct  "<< IsTr <<std::endl;
    
    mdl::Mesh new_res = {ver, indices, IsTr};
    ourModel.meshes[0] = new_res;


    // ................................................................................................................................................................

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // don't forget to enable shader before setting uniforms
        ourShader.use();

        // view/projection transformations
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);
        // hardcoded light
        ourShader.setVec3("objectColor", 1.0f, 1.f, 1.f);

        ourShader.setVec3("lightColor", 0.8f, 0.6f, 0.2f);
        ourShader.setVec3("lightDir", 1.f, 2.f, 3.f);

        ourShader.setVec3("lightColor1", 0.0f, 0.2f, 0.3f);
        ourShader.setVec3("lightDir1", -1.f, -2.f, -3.f);

        // render the loaded model
        glm::mat4 model = glm::mat4( 1.0f );
        model = glm::translate(model, glm::vec3(0.0f, 0.0f, 0.0f)); // translate it down so it's at the center of the scene
        model = glm::scale(model, glm::vec3(1.0f, 1.0f, 1.0f));	// it's a bit too big for our scene, so scale it down
        ourShader.setMat4("model", model);
        ourModel.Draw(ourShader);


        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

    

    return 0;
}








// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
        WireframeMode = !WireframeMode;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}