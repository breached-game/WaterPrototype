using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VelocityGrid : ScriptableObject
{
    private float[,,] cellsVelocities;
    private float[,,] newCellsVelocities;
    private Vector3Int[] basisVectors = new Vector3Int[2];

    public void Setup(int width, int height, int depth, float initVelocity, Vector3Int b)
    {
        cellsVelocities = new float[width, height, depth];
        newCellsVelocities = new float[width, height, depth];
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    cellsVelocities[x, y, z] = initVelocity;
                }
            }
        }
        basisVectors[1] = b;
        basisVectors[0] = -1 * b;
    }

    public void SetVelocity(int x, int y, int z, float v)
    {
        cellsVelocities[x, y, z] = v;
    }

    //CHANGED cells velocities to be the negative of the basis vector b, since the index is xi-1/2

    public void SetVelocities(Dictionary<Vector3Int, float> vs, Vector3Int pos)
    {
        cellsVelocities[pos.x, pos.y, pos.z] = vs[basisVectors[0]]; //was previously basisVectors[1]
        pos += basisVectors[1];
        cellsVelocities[pos.x, pos.y, pos.z] = vs[basisVectors[1]]; //was previously basisVectors[0]
    }

    public void SetNewVelocities(Dictionary<Vector3Int, float> vs, Vector3Int pos)
    {
        newCellsVelocities[pos.x, pos.y, pos.z] = vs[basisVectors[0]]; //was previously basisVector[1]
        pos += basisVectors[1];
        newCellsVelocities[pos.x, pos.y, pos.z] = vs[basisVectors[1]]; //was previously basisVector[0]
    }

    public void SetNewVelocity(float v, Vector3Int pos)
    {
        newCellsVelocities[pos.x, pos.y, pos.z] = v;
    }

    public Dictionary<Vector3Int, float> GetVelocities(Vector3Int pos)
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        velocities.Add(basisVectors[0], cellsVelocities[pos.x, pos.y, pos.z]);
        pos += basisVectors[1];
        velocities.Add(basisVectors[1], cellsVelocities[pos.x, pos.y, pos.z]);
        return velocities;
    }

    public Dictionary<Vector3Int, float> GetNewVelocities(Vector3Int pos)
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        velocities.Add(basisVectors[0], newCellsVelocities[pos.x, pos.y, pos.z]);
        pos += basisVectors[1];
        velocities.Add(basisVectors[1], newCellsVelocities[pos.x, pos.y, pos.z]);
        return velocities;
    }

    public float[,,] GetVelocityGrid()
    {
        return cellsVelocities;
    }

    public float[,,] GetNewVelocityGrid()
    {
        return newCellsVelocities;
    }

    public void UpdateVelocities()
    {
        cellsVelocities = newCellsVelocities;
    }
}
