using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public enum Contents
{
    Full,
    Surface,
    Empty
}

public class GridCell : MonoBehaviour
{
    private Vector3Int position;
    private Dictionary<Vector3Int,float> velocities;
    private Dictionary<Vector3Int, float> newVelocities;
    private float pressure;
    private float newPressure;
    private float divergence;
    private Contents cellContents;
    public GridCell(Vector3Int arg_position, float arg_pressure, Contents c)
    {
        position = arg_position;
        velocities = new Dictionary<Vector3Int, float>();
        pressure = arg_pressure;
        Debug.Log("Cell at " + arg_position + " created");
        divergence = 0f;
        cellContents = c;
    }

    public void SetVelocities(Dictionary<Vector3Int, float> arg_velocities)
    {
        velocities = arg_velocities;
    }

    public Dictionary<Vector3Int, float> GetVelocities()
    {
        return velocities;
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
        velocities = newVelocities;
    }

    public void SetNewVelocities(Dictionary<Vector3Int, float> arg_velocities)
    {
        newVelocities = arg_velocities;
    }

    public Dictionary<Vector3Int, float> GetNewVelocities()
    {
        return newVelocities;
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
