using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;

public class WaterGrid : MonoBehaviour
{
    private Grid water_grid;
    private GridColumn[,] gridArray;
    private int width = 60;
    private int height = 60;
    private int depth = 60;
    private float dz;
    private float dx;
    private float cellSize;
    public float dt = 0.05f;
    public GameObject columnMesh;
    public Vector3Int[] inflowLocations;

    private Vector3Int[] basisVectors = new Vector3Int[6] { Vector3Int.right, Vector3Int.left,
        Vector3Int.up, Vector3Int.down, Vector3Int.forward, Vector3Int.back};


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
                if (x > 5)
                {
                    column.Setup(new Vector2Int(x, z), 0, mesh);
                }
                else 
                { 
                    column.Setup(new Vector2Int(x, z), 0, mesh); 
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

        for (int i = 0; i < inflowLocationsSize; i++)
        {
            xInflow = inflowLocations[i][0];
            yInflow = inflowLocations[i][1];
            zInflow = inflowLocations[i][2];

            gridArray[xInflow, zInflow].Seth(yInflow);
        }
    }

    private float SumInflows(GridColumn currentColumn)
    {
        float iR;
        float iL;
        float iT;
        float iB;
        Vector2Int pos = currentColumn.GetPos();
        iL = gridArray[pos.x - 1, pos.y].GetNewOutflows()[Vector2Int.right];
        iR = gridArray[pos.x + 1, pos.y].GetNewOutflows()[Vector2Int.left];
        iT = gridArray[pos.x, pos.y + 1].GetNewOutflows()[Vector2Int.down];
        iB = gridArray[pos.x, pos.y - 1].GetNewOutflows()[Vector2Int.up];
        return iL + iR + iT + iB;
    }

    public void FixedUpdate()
    {
        float dhL, dhR, dhT, dhB;
        float dt_A_g_l = dt * Mathf.Pow(dx, 2) * 3f / dx;
        float K;
        float dV;
        float totalHeight;
        float totalFlux;
        Dictionary<Vector2Int, float> tempFlux = new Dictionary<Vector2Int, float>();
        GridColumn currentColumn;

        int xInflow;
        int yInflow;
        int zInflow;
        int inflowLocationsSize = inflowLocations.Length;

        for (int i = 0; i < inflowLocationsSize; i++)
        {
            xInflow = inflowLocations[i][0];
            yInflow = inflowLocations[i][1];
            zInflow = inflowLocations[i][2];

            gridArray[xInflow, zInflow].Seth(gridArray[xInflow, zInflow].Geth() + 100f * dt);
        }


        for (int x = 1; x < width - 1; x++)
        {
            for (int z = 1; z < height - 1; z++)
            {
                tempFlux.Clear();
                currentColumn = gridArray[x, z];
                totalHeight = currentColumn.GetH() + currentColumn.Geth();
                dhL = totalHeight - gridArray[x - 1, z].GetH() - gridArray[x-1, z].Geth();
                dhR = totalHeight - gridArray[x + 1, z].GetH() - gridArray[x + 1, z].Geth();
                dhT = totalHeight - gridArray[x, z + 1].GetH() - gridArray[x, z + 1].Geth();
                dhB = totalHeight - gridArray[x , z - 1].GetH() - gridArray[x, z -1].Geth();

                tempFlux.Add(Vector2Int.left, Mathf.Max(0.0f, currentColumn.GetOutflows()[Vector2Int.left] + dt_A_g_l * dhL));
                tempFlux.Add(Vector2Int.right, Mathf.Max(0.0f, currentColumn.GetOutflows()[Vector2Int.right] + dt_A_g_l * dhR));
                tempFlux.Add(Vector2Int.up, Mathf.Max(0.0f, currentColumn.GetOutflows()[Vector2Int.up] + dt_A_g_l * dhT));
                tempFlux.Add(Vector2Int.down, Mathf.Max(0.0f, currentColumn.GetOutflows()[Vector2Int.down] + dt_A_g_l * dhB));

                if (x == 1)
                {
                    tempFlux[Vector2Int.left] = 0f;
                }
                else if(x == width - 2)
                {
                    tempFlux[Vector2Int.right] = 0f;
                }
                if (z == 1)
                {
                    tempFlux[Vector2Int.down] = 0f;
                }
                else if(z == height - 2)
                {
                    tempFlux[Vector2Int.up] = 0f;
                }

                totalFlux = tempFlux.Sum(x => x.Value);
                if (totalFlux > 0.0f)
                {
                    K = Mathf.Min(1f, currentColumn.Geth() * dx * dx / totalFlux /dt);
                    tempFlux[Vector2Int.left] = K * tempFlux[Vector2Int.left];
                    tempFlux[Vector2Int.right] = K * tempFlux[Vector2Int.right];
                    tempFlux[Vector2Int.up] = K * tempFlux[Vector2Int.up];
                    tempFlux[Vector2Int.down] = K * tempFlux[Vector2Int.down];
                }
                currentColumn.SetNewOutflows(tempFlux);
            }
        }
        for (int x = 1; x < width - 1; x++)
        {
            for (int z = 1; z < height - 1; z++)
            {
                currentColumn = gridArray[x, z];
                dV = dt * (SumInflows(currentColumn) - currentColumn.GetNewOutflows().Sum(x => x.Value));
                // May want to change from int
                currentColumn.SetNewh(currentColumn.Geth() + dV / (dx * dx));            
            }
        }
        foreach (GridColumn column in gridArray)
        {
            column.UpdateValues();
            column.SetColumn();
        }

    }

}
