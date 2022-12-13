#version 330 core
out vec4 FragColor;

in vec3 Normal;
  
uniform vec3 lightDir; 
uniform vec3 lightColor;

uniform vec3 lightDir1; 
uniform vec3 lightColor1;

uniform vec3 objectColor;

void main()
{
    // ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;
  	
    // diffuse 
    vec3 ld = normalize(lightDir);
    vec3 norm = normalize(Normal);
    float diff = max(dot(norm, ld), 0.0);
    vec3 diffuse = diff * lightColor;

    vec3 ld1 = normalize(lightDir1);
    float diff1 = max(dot(norm, ld1), 0.0);
    vec3 diffuse1 = diff1 * lightColor1;
            
    vec3 result = (ambient + diffuse + diffuse1) * objectColor;
    FragColor = vec4(result, 1.0);
} 