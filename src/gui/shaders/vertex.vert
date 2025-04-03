#version 400 core
layout (location = 0) in vec3 aPos;
// layout (location = 1) in vec3 color; // Unused
layout (location = 2) in vec2 displacement;

out vec3 vertexColor; // Output color to fragment shader

// Scaling factor for the displacement effect on position
const float FACTOR = 1000.0f;

// --- TUNABLE PARAMETER ---
// Represents the scaled displacement magnitude that should result in full red.
// Green will appear at half this magnitude. Adjust as needed.
const float MAX_EXPECTED_SCALED_MAGNITUDE = 0.01;
// -------------------------

void main()
{
    // Calculate the displaced vertex position
    vec3 displacedPos = vec3(aPos.x + displacement.x * FACTOR, aPos.y + displacement.y * FACTOR, aPos.z);
    gl_Position = vec4(displacedPos, 1.0);

    // --- Color Calculation ---

    // Calculate the magnitude (length) of the scaled displacement vector
    float displacementMagnitude = length(displacement * FACTOR);

    // Normalize the magnitude to a [0.0, 1.0] range.
    float normalizedMagnitude = clamp(displacementMagnitude / MAX_EXPECTED_SCALED_MAGNITUDE, 0.0, 1.0);

    // Define the key colors for the gradient
    vec3 blueColor  = vec3(0.0, 0.0, 1.0);
    vec3 greenColor = vec3(0.0, 1.0, 0.0);
    vec3 redColor   = vec3(1.0, 0.0, 0.0);

    // Interpolate Blue -> Green -> Red
    if (normalizedMagnitude < 0.5)
    {
        // --- First half: Interpolate between Blue and Green ---
        // We need to map the range [0.0, 0.5) to [0.0, 1.0) for the mix factor.
        // factor = normalizedMagnitude / 0.5
        float factor = normalizedMagnitude * 2.0;
        vertexColor = mix(blueColor, greenColor, factor);
    }
    else
    {
        // --- Second half: Interpolate between Green and Red ---
        // We need to map the range [0.5, 1.0] to [0.0, 1.0] for the mix factor.
        // factor = (normalizedMagnitude - 0.5) / 0.5
        float factor = (normalizedMagnitude - 0.5) * 2.0;
        vertexColor = mix(greenColor, redColor, factor);
    }
}
