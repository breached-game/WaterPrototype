using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum Contents
{
    Full,
    Surface,
    Empty,
    Unknown,
    Solid
}

public class GridCell : ScriptableObject
{
    private Vector3Int position;
    private float pressure;
    private float newPressure;
    private float divergence;
    private Contents cellContents;
    private List<Particle> particles;
    public void Setup(Vector3Int arg_position, float arg_pressure, Contents c)
    {
        position = arg_position;
        pressure = arg_pressure;
        divergence = 0f;
        cellContents = c;
        particles = new List<Particle>();
    }

    public void AddParticle(Particle particle)
    {
        particles.Add(particle);
    }

    public void RemoveParticle(Particle particle)
    {
        particles.Remove(particle);
    }

    public List<Particle> GetParticles()
    {
        return particles;
    }

    public Vector3Int GetPos()
    {
        return position;
    }

    public float GetPressure()
    {
        return pressure;
    }

    public float GetDivergence()
    {
        return divergence;
    }
    public void SetDivergence(float d)
    {
        divergence = d;
    }

    public void UpdateValues()
    {
        pressure = newPressure;
    }

    public void SetNewPressure(float arg_presssure)
    {
        newPressure = arg_presssure;
    }

    public Contents GetContents()
    {
        return cellContents;
    }

    public void SetContents(Contents c)
    {
        cellContents = c;
    }
}
