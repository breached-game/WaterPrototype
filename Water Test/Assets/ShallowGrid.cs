using System.Collections;
using System.Collections.Generic;
using System;
using System.Linq;
using UnityEngine;

public class ShallowGrid : MonoBehaviour
{

    private VelocityGrid cellsVelocities_X;
    private VelocityGrid cellsVelocities_Z;


    private Grid water_grid;
    private GridColumn[,] gridArray;
    private int width = 10;
    private int height = 10;
    private int depth = 10;
    private float dz;
    private float dx;
    private float cellSize;
    private float gravity = 9.81f;
    public float dt = 0.05f;
    public Vector3Int[] inflowLocations;
    public float[] inflowVelocities;
    public Vector3Int[] inflowDirections;
    public GameObject columnMesh;

    private Vector3Int[] basisVectors = new Vector3Int[6] { Vector3Int.right, Vector3Int.left,
        Vector3Int.up, Vector3Int.down, Vector3Int.forward, Vector3Int.back};
    // Start is called before the first frame update
    void Awake()
    {
        GameObject mesh;
        Time.fixedDeltaTime = dt;
        water_grid = gameObject.GetComponent<Grid>();
        cellSize = water_grid.cellSize[0];
        dx = cellSize;
        dz = cellSize;
        gridArray = new GridColumn[width, depth];
        for (int x = 0; x < width; x++)
        {
            for (int z = 0; z < depth; z++)
            {
                GridColumn column = ScriptableObject.CreateInstance<GridColumn>();
                mesh = GameObject.Instantiate(columnMesh, water_grid.CellToLocal(new Vector3Int(x, 0, z)), Quaternion.identity);

                if (x == 0 || z == 0 || x == width - 1 || z == depth - 1)
                {
                    column.Setup(new Vector2Int(x, z), height, mesh);
                }
                else
                {
                    column.Setup(new Vector2Int(x, z), 1, mesh);
                }
                gridArray[x, z] = column;
            }
        }
        Setup();
    }

    private void Setup()
    {
        int xInflow;
        int yInflow;
        int zInflow;
        int inflowLocationsSize = inflowLocations.Length;
        int inflowDirectionsSize = inflowDirections.Length;
        int inflowVelocitiesSize = inflowVelocities.Length;

        cellsVelocities_X = ScriptableObject.CreateInstance<VelocityGrid>();
        cellsVelocities_Z = ScriptableObject.CreateInstance<VelocityGrid>();
        cellsVelocities_X.Setup(width + 1, depth, 0f, Vector2Int.right);
        cellsVelocities_Z.Setup(width, depth + 1, 0f, Vector2Int.up);

        if (!((inflowLocationsSize == inflowDirectionsSize) & (inflowLocationsSize == inflowVelocitiesSize)))
        {
            throw new Exception("Inflow sizes do not match");
        }

        for (int i = 0; i < inflowLocationsSize; i++)
        {
            xInflow = inflowLocations[i][0];
            yInflow = inflowLocations[i][1];
            zInflow = inflowLocations[i][2];

            // Need to add inflowDirection
            if (inflowDirections[i][0] != 0)
            {
                cellsVelocities_X.SetVelocity(xInflow, zInflow, (((inflowDirections[i][0] + 1) / 2) - 1) * inflowVelocities[i]);
                cellsVelocities_X.SetVelocity(xInflow + 1, zInflow, ((inflowDirections[i][0] + 1) / 2) * inflowVelocities[i]); // x+1/2
            }
            if (inflowDirections[i][1] != 0)
            {
                gridArray[xInflow, zInflow].Seth(yInflow - 1); //-1 cause H right now is 1
            }
            if (inflowDirections[i][2] != 0)
            {
                cellsVelocities_Z.SetVelocity(xInflow, zInflow, (((inflowDirections[i][2] + 1) / 2) - 1) * inflowVelocities[i]);
                cellsVelocities_Z.SetVelocity(xInflow, zInflow + 1, ((inflowDirections[i][2] + 1) / 2) * inflowVelocities[i]); // z+1/2
            }

            //gridArray[inflowLocations[i][0], inflowLocations[i][1], inflowLocations[i][2]].SetContents(Contents.Surface);
        }
    }

    private float EtaUpdate(GridColumn c)
    {
        float currentEta = c.Geth() + c.GetH();
        Vector2Int pos = c.GetPos();
        float uComponent = (cellsVelocities_X.GetVelocityGrid()[pos.x + 1, pos.y] - cellsVelocities_X.GetVelocityGrid()[pos.x, pos.y]) / dx;
        float wComponent = (cellsVelocities_Z.GetVelocityGrid()[pos.x, pos.y + 1] - cellsVelocities_Z.GetVelocityGrid()[pos.x, pos.y]) / dz;
        
        return -currentEta * (uComponent + wComponent);
    }

    private void SetVelocities()
    {
        float X_velocity;
        float Z_velocity;
        for (int x = 0; x < width - 1; x++)
        {
            for (int z = 0; z < depth - 1; z++)
            {
                X_velocity = -gravity * ((gridArray[x + 1, z].Geth() - gridArray[x, z].Geth()) / dx);
                Z_velocity = -gravity * ((gridArray[x, z].Geth() - gridArray[x, z + 1].Geth()) / dz);

                cellsVelocities_X.SetNewVelocity(X_velocity, new Vector2Int(x, z));
                cellsVelocities_Z.SetNewVelocity(Z_velocity, new Vector2Int(x, z));
            }
        }
    }

    public void FixedUpdate()
    {
        float newEta;
        float newh;
        SetVelocities();
        foreach (GridColumn c in gridArray)
        {
            newEta = EtaUpdate(c) * dt;
            newh = newEta - c.GetH();
            c.SetNewh((int) Math.Floor(newh));
            c.UpdateValues();
        }
        cellsVelocities_X.UpdateVelocities();
        cellsVelocities_Z.UpdateVelocities();
        foreach (GridColumn c in gridArray)
        {
            c.SetColumn();
        }
    }
}
