using System.Collections;
using System.Collections.Generic;
using System;
using UnityEngine;

public class WaterGrid : MonoBehaviour
{

    private Grid water_grid;
    private GridCell[,,] gridArray;
    private int width, height, depth = 10;
    private float dx, dy, dz = 0.2f;
    private float dt = 1f;
    private float vis = 1.3f;
    private float beta;
    private float b_0 = 1.7f;
    private float eps = 0.0001f;
    public Vector3 inflowLocation;
    public float inflowVelocity;
    public Vector3Int inflowDirection;
    [SerializeField] private GameObject WaterParticle;
    private Particle particle;

    private Vector3Int[] basisVectors = new Vector3Int[6] { new Vector3Int(1,0,0), new Vector3Int(-1,0,0), 
        new Vector3Int(0,1,0), new Vector3Int(0,-1,0), new Vector3Int(0,0,1), new Vector3Int(0,0,-1)};

    void Start()
    {
        Time.fixedDeltaTime = dt;
        water_grid = gameObject.GetComponent<Grid>();
        gridArray = new GridCell[width, height, depth];
        Setup();
        beta = (b_0 / (2 * dt)) * (1 / Mathf.Pow(dx, 2)) * (1 / Mathf.Pow(dy, 2)) * (1 / Mathf.Pow(dz, 2));
    }

    private void Setup()
    {
        float[] velocities = new float[6] {0f,0f,0f,0f,0f,0f};
        particle = new Particle(inflowLocation, WaterParticle);
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    GridCell cell = new GridCell(new Vector3Int(x, y, z), 0f, Contents.Full);
                    gridArray[x, y, z] = cell;
                    if (inflowLocation[0] == x & inflowLocation[1] == y & inflowLocation[2] == z)
                    {
                        velocities[Array.IndexOf(basisVectors, inflowDirection)] = inflowVelocity;
                    }
                    gridArray[x, y, z].SetVelocities(InstVelocities(velocities));
                }
            }
        }
    }

    private Dictionary<Vector3Int, float> InstVelocities(float[] vs)
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        for (int i = 0; i < 6; i++)
        {
            velocities.Add(basisVectors[i], vs[i]);
        }
        return velocities;
    }
    
    private Dictionary<Vector3Int,GridCell> GetNeighbours(GridCell c)
    {
        Vector3Int neighbour;
        Dictionary<Vector3Int, GridCell> neighbours = new Dictionary<Vector3Int, GridCell>();
        Vector3Int pos = c.GetPos();
        for (int i = 0; i < 6; i++)
        {
            neighbour = pos + basisVectors[i];
            if (!(neighbour[0] > width - 1 | neighbour[0] < 0))
            {
                if (!(neighbour[1] > height - 1 | neighbour[1] < 0))
                {
                    if (!(neighbour[2] > depth - 1 | neighbour[2] < 0))
                    {
                        neighbours.Add(basisVectors[i],gridArray[neighbour[0], neighbour[1], neighbour[2]]);
                    }
                }
            }
        }
        return neighbours;
    }

    private GridCell GetNeighbour(GridCell c, Dictionary<Vector3Int, GridCell> ns, Vector3Int d)
    {
        GridCell neighbour;
        if (ns.TryGetValue(d,out neighbour)) {
        }
        else
        {
            neighbour = new GridCell(new Vector3Int(-1, -1, -1), 0, Contents.Empty);
            neighbour.SetVelocities(InstVelocities(new float[] { 0f, 0f, 0f, 0f, 0f, 0f }));
        }
        return neighbour;
    }

    private float VelocityTilde(GridCell c, float v, Vector3Int d, Dictionary<Vector3Int, GridCell> ns)
    {
        Dictionary<Vector3Int, float> velocities;
        float neighbourVelocity;
        float neighbourAverage;
        float cellAverage;
        float c1;
        float g = 0f;
        float cx;
        float cy;
        float cz;
        float c2;
        float c3;
        float c4;
        float c5;
        float c6;

        velocities = c.GetVelocities();
        neighbourVelocity = GetNeighbour(c, ns, d).GetVelocities()[d];// u_i+1.5,j,k
        neighbourAverage = 0.5f * (neighbourVelocity + velocities[d]);// u_i+1,j,k
        cellAverage = 0.5f * (velocities[d] + velocities[-1 * d]);
        c1 = velocities[d] + dt * ((1 / dx) * (Mathf.Pow(cellAverage, 2) - Mathf.Pow(neighbourAverage, 2)));
        cx = (1 / dx) * ((velocities[d] * velocities[new Vector3Int(0, d[0] * -1, 0)]) - (velocities[d] * velocities[new Vector3Int(0, d[0], 0)]));
        cy = (1 / dy) * ((velocities[d] * velocities[new Vector3Int(d[0] * -1, 0, 0)]) - (velocities[d] * velocities[new Vector3Int(d[0], 0, 0)]));
        cz = (1 / dz) * ((velocities[d] * velocities[new Vector3Int(0, 0, d[0] * -1)]) - (velocities[d] * velocities[new Vector3Int(0, 0, d[0])]));
        c2 = 0f;
        c3 = 0f;
        if (d.Equals(basisVectors[0])) {
            c2 = cy + cz;
            c3 = (1 / dx) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[1])) {
            c2 = cy + cz;
            c3 = (1 / dx) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[2])) {
            c2 = cx + cz;
            g = 9.81f;
            c3 = (1 / dy) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[3])) {
            c2 = cx + cz;
            g = 9.81f;
            c3 = (1 / dy) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[4])) {
            c2 = cy + cx;
            c3 = (1 / dz) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[5])) {
            c2 = cy + cx;
            c3 = (1 / dz) * (c.GetPressure() - ns[d].GetPressure());
        }
        c4 = (vis / Mathf.Pow(dx, 2)) * (GetNeighbour(c, ns, basisVectors[0]).GetVelocities()[d] - 2 * v + GetNeighbour(c, ns, new Vector3Int(-1,0,0)).GetVelocities()[d]);
        c5 = (vis / Mathf.Pow(dy, 2)) * (GetNeighbour(c, ns, basisVectors[2]).GetVelocities()[d] - 2 * v + GetNeighbour(c, ns, basisVectors[3]).GetVelocities()[d]);
        c6 = (vis / Mathf.Pow(dz, 2)) * (GetNeighbour(c, ns, basisVectors[4]).GetVelocities()[d] - 2 * v + GetNeighbour(c, ns, basisVectors[5]).GetVelocities()[d]);
        return c1 + c2 + c3 + g + c4 + c5 + c6;
    }

    private float Divergence(GridCell c)
    {
        float cx;
        float cy;
        float cz;
        Dictionary<Vector3Int, float> velocities;
        velocities = c.GetVelocities();
        cx = (1 / dx) * (velocities[basisVectors[0]] - velocities[basisVectors[1]]);
        cy = (1 / dy) * (velocities[basisVectors[2]] - velocities[basisVectors[3]]);
        cz = (1 / dz) * (velocities[basisVectors[4]] - velocities[basisVectors[5]]);
        return -(cx + cy + cz);
    }

    private void VelocityUpdate(GridCell c)
    {
        float dp;
        Vector3Int b;
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = c.GetNewVelocities();
        float[] ds = new float[3] { dx, dy, dz };
        float d = Divergence(c);
        dp = beta * d;
        for (int i = 0; i < 6; i++)
        {
            b = basisVectors[i];
            velocities.Add(b,VelocityTilde(c, currentVelocities[b], b, GetNeighbours(c)));
            if (i % 2 == 0)
            {
                velocities[b] = velocities[b] + ((dt / (ds[Decimal.ToInt64(Decimal.Truncate(i / 2))])) * dp);
            }
            else
            {
                velocities[b] = velocities[b] - ((dt / (ds[Decimal.ToInt64(Decimal.Truncate(i / 2))])) * dp);
            }
        }
        c.SetDivergence(d);
        c.SetNewVelocities(velocities);
        c.SetNewPressure(c.GetPressure() + dp);
    }

    private void PressureIterations()
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Vector3Int b;
        bool con = true;
        // Runs VelocityUpdate after previously calculating Velocity Tilde for each cell
        // Keeps running it until divergence is less than epsilon (pre-defined)
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (gridArray[x, y, z].GetContents() == Contents.Full)
                    {
                        velocities.Clear();
                        for (int i = 0; i < 6; i++)
                        {
                            b = basisVectors[i];
                            velocities.Add(b, VelocityTilde(gridArray[x, y, z], gridArray[x, y, z].GetVelocities()[b], b, GetNeighbours(gridArray[x, y, z])));
                        }
                        gridArray[x, y, z].SetNewVelocities(velocities);
                    }
                }
            }
        }
        while (con)
        {
            con = false;
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    for (int z = 0; z < depth; z++)
                    {
                        VelocityUpdate(gridArray[x, y, z]);
                        if (gridArray[x, y, z].GetDivergence() > eps)
                        {
                            con = true;
                        }
                    }
                }
            }
        }
        Debug.Log("Cycle Finished");
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    gridArray[x, y, z].UpdateValues();
                }
            }
        }
    }

    private void ParticleVelocityUpdate(Particle p)
    {
        Vector3 position = p.getPosition();
        Vector3Int cellPosition = water_grid.LocalToCell(position);
        GridCell cell = gridArray[cellPosition[0], cellPosition[1], cellPosition[2]];


    }

    private void ParticleLocationUpdate(Particle p, Vector3 v)
    {
        Vector3 newPos = p.getPosition() + (v * dt);
        p.UpdatePosition(newPos);
    }

    public void FixedUpdate()
    {
        PressureIterations();
    }
}
