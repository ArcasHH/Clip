#version 330 core
out vec4 FragColor;

in vec3 Normal;
  
uniform vec3 lightDir; 
uniform vec3 lightColor;
uniform vec3 objectColor;

void main()
{
    // ambient
    float ambientStrength = 0.45;
    vec3 ambient = ambientStrength * lightColor;
  	
    // diffuse 
    vec3 ld = normalize(lightDir);
    vec3 norm = normalize(Normal);
    float diff = max(dot(norm, ld), 0.0);
    vec3 diffuse = diff * lightColor;
            
    vec3 result = (ambient + diffuse) * objectColor;
    FragColor = vec4(result, 1.0);
} 