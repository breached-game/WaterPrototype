using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particle : MonoBehaviour
{
    private Vector3Int position;
    private GameObject WaterParticle;
    // Start is called before the first frame update
    public Particle(Vector3Int initialPosition, GameObject p)
    {
        position = initialPosition;
        WaterParticle = p;
        WaterParticle.transform.position = position;
    }

    // Update is called once per frame
    void UpdatePosition(Vector3Int newPos)
    {
        position = newPos;
        WaterParticle.transform.position = position;
    }

    public Vector3Int getPosition()
    {
        return position;
    }
}
