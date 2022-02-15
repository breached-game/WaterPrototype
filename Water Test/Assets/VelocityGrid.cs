using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VelocityGrid : ScriptableObject
{
    private float[,] cellsVelocities;
    private float[,] newCellsVelocities;
    private Vector2Int[] basisVectors = new Vector2Int[2];

    public void Setup(int width, int height, float initVelocity, Vector2Int b)
    {
        cellsVelocities = new float[width, height];
        newCellsVelocities = new float[width, height];
        for (int x = 0; x < width; x++)
        {
            for (int z = 0; z < height; z++)
            {
                cellsVelocities[x, z] = initVelocity;
            }
            
        }
        basisVectors[1] = b;
        basisVectors[0] = -1 * b;
    }

    public void SetVelocity(int x, int z, float v)
    {
        cellsVelocities[x, z] = v;
    }

    //CHANGED cells velocities to be the negative of the basis vector b, since the index is xi-1/2

    public void SetVelocities(Dictionary<Vector2Int, float> vs, Vector2Int pos)
    {
        cellsVelocities[pos.x, pos.y] = vs[basisVectors[0]]; //was previously basisVectors[1]
        pos += basisVectors[1];
        cellsVelocities[pos.x, pos.y] = vs[basisVectors[1]]; //was previously basisVectors[0]
    }

    public void SetNewVelocities(Dictionary<Vector2Int, float> vs, Vector2Int pos)
    {
        newCellsVelocities[pos.x, pos.y] = vs[basisVectors[0]]; //was previously basisVector[1]
        pos += basisVectors[1];
        newCellsVelocities[pos.x, pos.y] = vs[basisVectors[1]]; //was previously basisVector[0]
    }

    public void SetNewVelocity(float v, Vector2Int pos)
    {
        newCellsVelocities[pos.x, pos.y] = v;
    }

    public Dictionary<Vector2Int, float> GetVelocities(Vector2Int pos)
    {
        Dictionary<Vector2Int, float> velocities = new Dictionary<Vector2Int, float>();
        velocities.Add(basisVectors[0], cellsVelocities[pos.x, pos.y]);
        pos += basisVectors[1];
        velocities.Add(basisVectors[1], cellsVelocities[pos.x, pos.y]);
        return velocities;
    }

    public Dictionary<Vector2Int, float> GetNewVelocities(Vector2Int pos)
    {
        Dictionary<Vector2Int, float> velocities = new Dictionary<Vector2Int, float>();
        velocities.Add(basisVectors[0], newCellsVelocities[pos.x, pos.y]);
        pos += basisVectors[1];
        velocities.Add(basisVectors[1], newCellsVelocities[pos.x, pos.y]);
        return velocities;
    }

    public float[,] GetVelocityGrid()
    {
        return cellsVelocities;
    }

    public float[,] GetNewVelocityGrid()
    {
        return newCellsVelocities;
    }

    public void UpdateVelocities()
    {
        cellsVelocities = newCellsVelocities;
    }
}
