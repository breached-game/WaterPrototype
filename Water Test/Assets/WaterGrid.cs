using System.Collections;
using System.Collections.Generic;
using System;
using System.Linq;
using UnityEngine;

public class WaterGrid : MonoBehaviour
{

    struct VelocityDistance
    {
        public float[] distance;
        public float[] velocity;
    }

    private Grid water_grid;
    private GridCell[,,] gridArray;
    private int width = 10;
    private int height = 10;
    private int depth = 10;
    private float dz = 0.2f;
    private float dx = 0.2f;
    private float dy = 0.2f;
    private float cellSize = 0.2f;
    private float dt = 0.1f;
    private float vis = 1.3f;
    private float beta;
    private float b_0 = 1.7f;
    private float eps = 0.0001f;
    public Vector3Int inflowLocation;
    public float inflowVelocity;
    public Vector3Int inflowDirection;
    [SerializeField] private GameObject WaterParticle;
    private Particle particle;

    private Vector3Int[] basisVectors = new Vector3Int[6] { new Vector3Int(1,0,0), new Vector3Int(-1,0,0), 
        new Vector3Int(0,1,0), new Vector3Int(0,-1,0), new Vector3Int(0,0,1), new Vector3Int(0,0,-1)};

    void Awake()
    {
        Time.fixedDeltaTime = dt;
        water_grid = gameObject.GetComponent<Grid>();
        gridArray = new GridCell[width, height, depth];
        Setup();
        beta = (b_0 / (2 * dt)) * (1 / Mathf.Pow(dx, 2)) * (1 / Mathf.Pow(dy, 2)) * (1 / Mathf.Pow(dz, 2));
        Debug.Log("Awake");
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (gridArray[x, y, z].GetVelocities() == null)
                    {
                        Debug.Log(gridArray[x, y, z].GetPos());
                    }
                }
            }
        }
    }

    void Start()
    {
        Debug.Log("Started");
    }

    private void Setup()
    {
        float[] velocities = new float[6] {0f,0f,0f,0f,0f,0f};
        particle = ScriptableObject.CreateInstance<Particle>();
        particle.Setup(water_grid.CellToLocal(inflowLocation), WaterParticle);
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    GridCell cell = ScriptableObject.CreateInstance<GridCell>();
                    cell.Setup(new Vector3Int(x, y, z), 0f, Contents.Empty);
                    gridArray[x, y, z] = cell;
                    if (inflowLocation[0] == x & inflowLocation[1] == y & inflowLocation[2] == z)
                    {
                        velocities[Array.IndexOf(basisVectors, inflowDirection)] = inflowVelocity;
                    }
                    gridArray[x, y, z].SetVelocities(InstVelocities(velocities));
                    if (x == 0 || y == 0 || z == 0 || x == width -1 || y == height-1 || z == depth -1)
                    {
                        gridArray[x, y, z].SetContents(Contents.Solid);
                    }
                }
            }
        }
        gridArray[inflowLocation[0], inflowLocation[1], inflowLocation[2]].SetContents(Contents.Surface);
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
            neighbour = ScriptableObject.CreateInstance<GridCell>();
            neighbour.Setup(new Vector3Int(-1, -1, -1), 0, Contents.Empty);
            neighbour.SetVelocities(InstVelocities(new float[] { 0f, 0f, 0f, 0f, 0f, 0f }));
        }
        for (int i = 0; i < width; i++)
        {
            for (int j = 0; j < height; j++)
            {
                for (int k = 0; k < depth; k++)
                {
                    if (gridArray[i, j, k].GetVelocities() == null)
                    {
                        Debug.Log(gridArray[i, j, k].GetPos());
                    }
                }
            }
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
        int sign = 1;

        Debug.Log("Velocity tilde called");

        velocities = c.GetVelocities();
        neighbourVelocity = GetNeighbour(c, ns, d).GetVelocities()[d];// u_i+1.5,j,k
        neighbourAverage = 0.5f * (neighbourVelocity + velocities[d]);// u_i+1,j,k
        cellAverage = 0.5f * (velocities[d] + velocities[-1 * d]);
        c1 = velocities[d] + (dt * ((1 / dx) * (Mathf.Pow(cellAverage, 2) - Mathf.Pow(neighbourAverage, 2))));
        if (d[0] + d[1] + d[2] < 0)
        {
            sign = -1;
        }
        cx = (1 / dx) * ((velocities[d] * velocities[new Vector3Int(0, sign * -1, 0)]) - (velocities[d] * velocities[new Vector3Int(0, sign, 0)]));
        cy = (1 / dy) * ((velocities[d] * velocities[new Vector3Int(sign * -1, 0, 0)]) - (velocities[d] * velocities[new Vector3Int(sign, 0, 0)]));
        cz = (1 / dz) * ((velocities[d] * velocities[new Vector3Int(0, 0, sign * -1)]) - (velocities[d] * velocities[new Vector3Int(0, 0, sign)]));
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
            g = -9.81f;
            c3 = (1 / dy) * (c.GetPressure() - ns[d].GetPressure());
        }
        else if (d.Equals(basisVectors[3])) {
            c2 = cx + cz;
            g = -9.81f;
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
        Debug.Log(c1 + c2 + c3 + g + c4 + c5 + c6);
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
            if (i % 2 == 0)
            {
                velocities[b] = currentVelocities[b] + ((dt / (ds[Decimal.ToInt64(Decimal.Truncate(i / 2))])) * dp);
            }
            else
            {
                velocities[b] = currentVelocities[b] - ((dt / (ds[Decimal.ToInt64(Decimal.Truncate(i / 2))])) * dp);
            }
        }
        c.SetDivergence(d);
        c.SetNewVelocities(velocities);
        c.SetNewPressure(c.GetPressure() + dp);
    }

    private void SurfaceVelocityUpdate(GridCell c)
    {
        float p = 0f;
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = c.GetVelocities();
        Dictionary<Vector3Int,GridCell> neighbours = GetNeighbours(c);
        GridCell neighbour;
        GridCell oppositeNeighbour;
        Contents contents;
        for (int i = 0; i < 6; i++)
        {
            Debug.Log("Cell: " + c.GetPos());
            Debug.Log(basisVectors[i]);
            neighbour = neighbours[basisVectors[i]];
            oppositeNeighbour = neighbours[basisVectors[i] * -1];
            contents = neighbour.GetContents();
            if ((contents.Equals(Contents.Full) || contents.Equals(Contents.Surface)) && !(oppositeNeighbour.GetContents().Equals(Contents.Full) || oppositeNeighbour.GetContents().Equals(Contents.Surface)))
            {
                velocities.Add(basisVectors[i] * -1, neighbour.GetVelocities()[basisVectors[i] * -1]);
            }
            else if (contents.Equals(Contents.Solid)) {
                velocities.Add(basisVectors[i], 0f);
            }
            else
            {
                velocities.Add(basisVectors[i], currentVelocities[basisVectors[i]]);
            }
        }
        c.SetNewVelocities(velocities);
        c.SetNewPressure(p);
    }

    private void SolidVelocityUpdate(GridCell c)
    {
        float p = 0f;
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = c.GetVelocities();
        Dictionary<Vector3Int, GridCell> neighbours = GetNeighbours(c);
        Dictionary<Vector3Int, float> neighbourVelocities;
        GridCell neighbour;

        for (int i = 0; i < 6; i++)
        {
            if (neighbours.TryGetValue(basisVectors[i], out neighbour))
            {
                if (!neighbour.GetContents().Equals(Contents.Solid))
                {
                    neighbourVelocities = neighbour.GetVelocities();
                    p = neighbour.GetPressure();
                    for (int j = 0; j < 6; j++)
                    { 
                        if (basisVectors[j] != (basisVectors[i] * -1))
                        {
                            velocities.Add(basisVectors[j], neighbourVelocities[basisVectors[j]]);
                        }
                        else
                        {
                            velocities.Add(basisVectors[j], 0f);
                        }
                    }
                    break;
                }
            }
            else
            {
            }
        }
        if (velocities.Count == 0)
        {
            velocities = InstVelocities(new float[] { 0f, 0f, 0f, 0f, 0f, 0f });
        }
        c.SetNewVelocities(velocities);
        c.SetNewPressure(p);
    }

    private void PressureIterations()
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Vector3Int b;
        bool con = true;
        // Runs VelocityUpdate after previously calculating Velocity Tilde for each cell
        // Keeps running it until divergence is less than epsilon (pre-defined)
        Debug.Log("Cycle Started");
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
        Debug.Log("Full cells set");
        while (con)
        {
            con = false;
            for (int x = 0; x < width; x++)
            {
                for (int y = 0; y < height; y++)
                {
                    for (int z = 0; z < depth; z++)
                    {
                        if (gridArray[x, y, z].GetContents().Equals(Contents.Full))
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
        }
        Debug.Log("Cycle Finished");
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (gridArray[x, y, z].GetContents() == Contents.Surface)
                    {
                        SurfaceVelocityUpdate(gridArray[x, y, z]);
                    }
                    else if (gridArray[x, y, z].GetContents() == Contents.Solid)
                    {
                        SolidVelocityUpdate(gridArray[x, y, z]);
                    }
                    else if (gridArray[x, y, z].GetContents() == Contents.Empty)
                    {
                        gridArray[x, y, z].SetNewVelocities(gridArray[x,y,z].GetVelocities());
                    }
                    gridArray[x, y, z].UpdateValues();
                }
            }
        }
    }

    private VelocityDistance GetParticleVelocityComponents(Vector3 particlePosition, GridCell cell, Vector3 localCellPosition, Vector3 v1,Vector3Int v2,Vector3 v3,Vector3Int v4,Vector3 v5,Vector3Int v6)
    {
        Vector3[] velocityPositions = new Vector3[4];
        VelocityDistance componentVelocities = new VelocityDistance();
        componentVelocities.distance = new float[4];
        componentVelocities.velocity = new float[4];
        velocityPositions[0] = localCellPosition + (v1 * cellSize / 2);
        componentVelocities.velocity[0] = cell.GetVelocities()[v2];
        velocityPositions[1] = localCellPosition + (v3 * cellSize / 2);
        componentVelocities.velocity[1] = cell.GetVelocities()[v4];
        velocityPositions[2] = localCellPosition + (v5 * cellSize) + (v1 * cellSize / 2);
        componentVelocities.velocity[2] = GetNeighbour(cell, GetNeighbours(cell), v6).GetVelocities()[v2];
        velocityPositions[3] = localCellPosition + (v5 * cellSize) + (v3 * cellSize / 2);
        componentVelocities.velocity[3] = GetNeighbour(cell, GetNeighbours(cell), v6).GetVelocities()[v4];
        for (int i = 0; i < 4; i++)
        {
            componentVelocities.distance[i] = Vector3.Magnitude(particlePosition - velocityPositions[i]);
        }
        return componentVelocities;
    } 

    private Vector3 ParticleVelocityUpdate(Particle p)
    {
        // Calculated in 2D so x + y, y + x, z + y
        Vector3 position = p.getPosition();
        float[] newVelocities = new float[3];
        int index;
        VelocityDistance[] componentsVelocities = new VelocityDistance[3];
        float distance;
        Vector3Int cellPosition = water_grid.LocalToCell(position);
        Debug.Log(cellPosition);
        GridCell cell = gridArray[cellPosition[0], cellPosition[1], cellPosition[2]];
        Vector3 localCellPosition = water_grid.CellToLocal(cellPosition) + Vector3.one * dx/2;
        //x component
        index = 1;
        distance = (position - localCellPosition)[index];
        if (distance < cellSize / 2)
        {
            componentsVelocities[0] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.left, Vector3Int.left, Vector3.right, Vector3Int.right, Vector3.down, Vector3Int.down);     
        }
        else
        {
            componentsVelocities[0] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.left, Vector3Int.left, Vector3.right, Vector3Int.right, Vector3.up, Vector3Int.up);
        }
        //y component
        index = 0;
        distance = (position - localCellPosition)[index];
        if (distance < cellSize / 2)
        {
            componentsVelocities[1] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.up, Vector3Int.up, Vector3.down, Vector3Int.down, Vector3.left, Vector3Int.left);
        }
        else
        {
            componentsVelocities[1] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.up, Vector3Int.up, Vector3.down, Vector3Int.down, Vector3.right, Vector3Int.right);

        }
        //z component
        index = 0;
        distance = (position - localCellPosition)[index];
        if (distance < cellSize / 2)
        {
            componentsVelocities[2] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.back, Vector3Int.back, Vector3.forward, Vector3Int.forward, Vector3.left, Vector3Int.left);
        }
        else
        {
            componentsVelocities[2] = GetParticleVelocityComponents(position, cell, localCellPosition, Vector3.back, Vector3Int.back, Vector3.forward, Vector3Int.forward, Vector3.right, Vector3Int.right);
        }
        for (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 4; i++)
            {
                newVelocities[j] += componentsVelocities[j].velocity[i] * componentsVelocities[j].distance[i];
            }
            newVelocities[j] = newVelocities[j] / (componentsVelocities[j].distance.Sum());
        }
        return new Vector3(newVelocities[0], newVelocities[1], newVelocities[2]);
    }

    private void ParticleLocationUpdate(Particle p)
    {
        Vector3 v;
        Vector3 newPos;
        Vector3Int cellPos;
        Vector3Int particleCellPos = water_grid.LocalToCell(p.getPosition());
        v = ParticleVelocityUpdate(p);
        Debug.Log("Velocity: " + v);
        newPos = p.getPosition() + (v * dt);
        cellPos = water_grid.LocalToCell(newPos);
        if (cellPos[0] == height-1 || cellPos[1] == width-1 || cellPos[2] == depth - 1)
        {
            newPos = p.getPosition();
            cellPos = water_grid.LocalToCell(newPos);
        }
        Debug.Log("cell pos: " + cellPos);
        gridArray[particleCellPos.x, particleCellPos.y, particleCellPos.z].RemoveParticle(p);
        gridArray[cellPos.x, cellPos.y, cellPos.z].AddParticle(p);
        p.UpdatePosition(newPos);
    }

    public void CellContents()
    {
        Dictionary<Vector3Int, GridCell> neighbours;
        // MAKE MORE EFFICIENT
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (!gridArray[x, y, z].GetContents().Equals(Contents.Solid) & gridArray[x,y,z].GetParticles().Count == 0) {
                        gridArray[x, y, z].SetContents(Contents.Empty);
                    }
                }
            }
        }
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (gridArray[x,y,z].GetParticles().Count > 0)
                    {

                        neighbours = GetNeighbours(gridArray[x,y,z]);
                        foreach (KeyValuePair<Vector3Int, GridCell> entry in neighbours)
                        {
                            if (entry.Value.GetContents().Equals(Contents.Empty))
                            {
                                gridArray[x, y, z].SetContents(Contents.Surface);
                                break;
                            }
                        }

                    }
                }
            }
        }
        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    if (gridArray[x,y,z].GetParticles().Count > 0 & !gridArray[x,y,z].GetContents().Equals(Contents.Surface) & !gridArray[x, y, z].GetContents().Equals(Contents.Solid))
                    {
                        gridArray[x, y, z].SetContents(Contents.Full);
                    }
                }
            }
        }
    }

    public void FixedUpdate()
    {
        CellContents();
        Debug.Log("Update");
        PressureIterations();
        ParticleLocationUpdate(particle);
    }
}
