using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GridCell : MonoBehaviour
{
    private Vector3Int position;
    private Dictionary<Vector3Int,float> velocities;
    private Dictionary<Vector3Int, float> newVelocities;
    private float pressure;
    private float newPressure;
    private float divergence;
    public GridCell(Vector3Int arg_position, float arg_pressure)
    {
        position = arg_position;
        velocities = new Dictionary<Vector3Int, float>();
        pressure = arg_pressure;
        Debug.Log("Cell at " + arg_position + " created");
        divergence = 0f;
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

    public void SetNewPressure(float arg_presssure)
    {
        newPressure = arg_presssure;
    }
}
