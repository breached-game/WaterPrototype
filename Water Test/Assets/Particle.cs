using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particle : MonoBehaviour
{
    private Vector3 position;
    private GameObject WaterParticle;
    // Start is called before the first frame update
    public Particle(Vector3 initialPosition, GameObject p)
    {
        position = initialPosition;
        WaterParticle = p;
        WaterParticle.transform.position = position;
    }

    // Update is called once per frame
    void UpdatePosition(Vector3 newPos)
    {
        position = newPos;
        WaterParticle.transform.position = position;
    }

    public Vector3 getPosition()
    {
        return position;
    }
}
