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

    private VelocityGrid cellsVelocities_X;
    private VelocityGrid cellsVelocities_Y;
    private VelocityGrid cellsVelocities_Z;


    private Grid water_grid;
    private GridCell[,,] gridArray;
    private int width = 10;
    private int height = 10;
    private int depth = 10;
    private float dz;
    private float dx;
    private float dy;
    private float cellSize;
    public float dt = 0.05f;
    public float vis = 1.3f;
    private float beta;
    private float b_0 = 1.7f;
    private float eps = 0.0001f;
    public Vector3Int inflowLocation;
    public float inflowVelocity;
    public Vector3Int inflowDirection;
    [SerializeField] private GameObject WaterParticle;
    private List<Particle> particleList = new List<Particle>();
    private int particleCount;

    private Vector3Int[] basisVectors = new Vector3Int[6] { new Vector3Int(1,0,0), new Vector3Int(-1,0,0),
        new Vector3Int(0,1,0), new Vector3Int(0,-1,0), new Vector3Int(0,0,1), new Vector3Int(0,0,-1)};

    void Awake()
    {
        Time.fixedDeltaTime = dt;
        water_grid = gameObject.GetComponent<Grid>();
        cellSize = water_grid.cellSize[0];
        dx = cellSize;
        dy = cellSize;
        dz = cellSize;
        gridArray = new GridCell[width, height, depth];
        Setup();
        beta = (b_0 / (2 * dt)) * (1 / Mathf.Pow(dx, 2)) * (1 / Mathf.Pow(dy, 2)) * (1 / Mathf.Pow(dz, 2));
    }

    private void Setup()
    {
        float gravityVelocity = (float)(-9.81 * dt);

        for (int x = 0; x < width; x++)
        {
            for (int y = 0; y < height; y++)
            {
                for (int z = 0; z < depth; z++)
                {
                    GridCell cell = ScriptableObject.CreateInstance<GridCell>(); //
                    cell.Setup(new Vector3Int(x, y, z), 0f, Contents.Empty);    //
                    gridArray[x, y, z] = cell;

                    if (x == 0 || y == 0 || z == 0 || x == width - 1 || y == height - 1 || z == depth - 1)
                    {
                        gridArray[x, y, z].SetContents(Contents.Solid);
                    }
                }
            }
        }

        cellsVelocities_X = ScriptableObject.CreateInstance<VelocityGrid>();
        cellsVelocities_Y = ScriptableObject.CreateInstance<VelocityGrid>();
        cellsVelocities_Z = ScriptableObject.CreateInstance<VelocityGrid>();
        cellsVelocities_X.Setup(width + 1, height, depth, 0f, Vector3Int.right);
        cellsVelocities_Y.Setup(width, height + 1, depth, gravityVelocity, Vector3Int.up);
        cellsVelocities_Z.Setup(width, height, depth + 1, 0f, Vector3Int.forward);
        int xInflow = inflowLocation[0];
        int yInflow = inflowLocation[1];
        int zInflow = inflowLocation[2];

        // Need to add inflowDirection
        cellsVelocities_X.SetVelocity(xInflow, yInflow, zInflow, (((inflowDirection[0] + 1) / 2) - 1) * inflowVelocity);
        cellsVelocities_X.SetVelocity(xInflow + 1, yInflow, zInflow, ((inflowDirection[0] + 1)/2) * inflowVelocity); // x+1/2
        cellsVelocities_Y.SetVelocity(xInflow, yInflow, zInflow, (((inflowDirection[1] + 1) / 2) - 1) * inflowVelocity);
        cellsVelocities_Y.SetVelocity(xInflow, yInflow + 1, zInflow, ((inflowDirection[1] + 1) / 2) * inflowVelocity); // x+1/2
        cellsVelocities_Z.SetVelocity(xInflow, yInflow, zInflow, (((inflowDirection[2] + 1) / 2) - 1) * inflowVelocity);
        cellsVelocities_Z.SetVelocity(xInflow, yInflow, zInflow + 1, ((inflowDirection[2] + 1) / 2) * inflowVelocity); // x+1/2

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

    private void CreateParticle(Vector3Int position)
    {
        Particle p;
        Vector3 localPosition;
        localPosition = water_grid.CellToLocal(position);
        GameObject particleObject = Instantiate(WaterParticle, localPosition, Quaternion.identity);
        p = ScriptableObject.CreateInstance<Particle>();
        p.Setup(water_grid.CellToLocal(position), particleObject);
        particleList.Add(p);
    }
    private Dictionary<Vector3Int, GridCell> GetNeighbours(GridCell cell)
    {
        // Passes in cell coordinates (could be changed to Vector3Int of position), returns dict of form -> direction from cell : position of neighbour

        GridCell neighbour;
        Vector3Int b; // Basis vector
        Vector3Int pos = cell.GetPos();
        Vector3Int neighbourPos;

        Dictionary<Vector3Int, GridCell> neighbours = new Dictionary<Vector3Int, GridCell>();


        for (int i = 0; i < 6; i++)
        {
            b = basisVectors[i];
            neighbourPos = pos + b;

            if (!(neighbourPos.x > width - 1 | neighbourPos.x < 0))
            {
                if (!(neighbourPos.y > height - 1 | neighbourPos.y < 0))
                {
                    if (!(neighbourPos.z > depth - 1 | neighbourPos.z < 0))
                    {
                        neighbour = gridArray[pos.x + b.x, pos.y + b.y, pos.z + b.z];
                        neighbours.Add(b, neighbour);
                    }
                }
            }
        }
        return neighbours;
    }

    private Dictionary<Vector3Int, float> GetCellAllVelocities(GridCell c)
    {
        Vector3Int pos = c.GetPos();
        Dictionary<Vector3Int, float> allVelocities = new Dictionary<Vector3Int, float>();
        cellsVelocities_X.GetVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        cellsVelocities_Y.GetVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        cellsVelocities_Z.GetVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        return allVelocities;
    }

    private Dictionary<Vector3Int, float> GetCellAllNewVelocities(GridCell c)
    {
        Vector3Int pos = c.GetPos();
        Dictionary<Vector3Int, float> allVelocities = new Dictionary<Vector3Int, float>();
        cellsVelocities_X.GetNewVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        cellsVelocities_Y.GetNewVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        cellsVelocities_Z.GetNewVelocities(c.GetPos()).ToList().ForEach(x => allVelocities.Add(x.Key, x.Value));
        return allVelocities;
    }



    private float GetNeighbourVelocity(Dictionary<Vector3Int, GridCell> neighbours, Vector3Int d)
    {
        //Note, will get the neighbour in direction d, and returns the neighbours velocity IN THE DIRECTION OF d !
        //Haven't developed a way of getting neighbour in direction d, and then selecting a velocity of that neighbour not in direction d.

        Vector3Int neighbour = neighbours[d].GetPos();
        return GetCellAllVelocities(gridArray[neighbour.x, neighbour.y, neighbour.z])[d]; //Could be done more efficiently,  would need a quick method for selecting which velocity array (either x, y, z) to choose from using just d
    }
    private float VelocityTilde(GridCell c, float v, Vector3Int d, Dictionary<Vector3Int, GridCell> neighbours)
    {

        float selectedVelocity;
        float pressure;
        float neighbourPressure;
        Vector3Int pos = c.GetPos();

        Dictionary<Vector3Int, float> allCellVelocities;
        Vector3Int neighbour;
        float neighbourVelocity;
        float neighbourAverage;
        float cellAverage;
        float g = 9.81f;
        float cx;
        float cy;
        float cz;
        float c1;
        float c2;
        float c3;
        float c4;
        float c5;
        float c6;
        int sign = 1;

        Debug.Log("Velocity tilde called");

        //Current cell
        allCellVelocities = GetCellAllVelocities(gridArray[pos.x, pos.y, pos.z]);
        selectedVelocity = allCellVelocities[d];
        pressure = gridArray[pos.x, pos.y, pos.z].GetPressure();

        //Neighbour cell
        neighbour = neighbours[d].GetPos();
        neighbourVelocity = GetCellAllVelocities(gridArray[neighbour.x, neighbour.y, neighbour.z])[d];
        neighbourPressure = gridArray[neighbour.x, neighbour.y, neighbour.z].GetPressure();


        neighbourAverage = 0.5f * (selectedVelocity + neighbourVelocity);
        cellAverage = 0.5f * (selectedVelocity + allCellVelocities[d * -1]);
        c1 = selectedVelocity + (dt * ((1 / dx) * (Mathf.Pow(cellAverage, 2) - Mathf.Pow(neighbourAverage, 2))));

        if (d[0] + d[1] + d[2] < 0)
        {
            sign = -1;
        }
        cx = (1 / dx) * ((selectedVelocity * allCellVelocities[new Vector3Int(0, sign * -1, 0)]) - (selectedVelocity * allCellVelocities[new Vector3Int(0, sign, 0)]));
        cy = (1 / dy) * ((selectedVelocity * allCellVelocities[new Vector3Int(sign * -1, 0, 0)]) - (selectedVelocity * allCellVelocities[new Vector3Int(sign, 0, 0)]));
        cz = (1 / dz) * ((selectedVelocity * allCellVelocities[new Vector3Int(0, 0, sign * -1)]) - (selectedVelocity * allCellVelocities[new Vector3Int(0, 0, sign)]));
        c2 = 0f;
        c3 = 0f;
        if (d.Equals(basisVectors[0]))
        {
            c2 = cy + cz;
            c3 = (1 / dx) * (pressure - neighbourPressure);
        }
        else if (d.Equals(basisVectors[1]))
        {
            c2 = cy + cz;
            c3 = (1 / dx) * (pressure - neighbourPressure);
        }
        else if (d.Equals(basisVectors[2]))
        {
            c2 = cx + cz;
            g = -9.81f;
            c3 = (1 / dy) * (pressure - neighbourPressure);
        }
        else if (d.Equals(basisVectors[3]))
        {
            c2 = cx + cz;
            g = -9.81f;
            c3 = (1 / dy) * (pressure - neighbourPressure);
        }
        else if (d.Equals(basisVectors[4]))
        {
            c2 = cy + cx;
            c3 = (1 / dz) * (pressure - neighbourPressure);
        }
        else if (d.Equals(basisVectors[5]))
        {
            c2 = cy + cx;
            c3 = (1 / dz) * (pressure - neighbourPressure);
        }




        c4 = (vis / Mathf.Pow(dx, 2)) * GetNeighbourVelocity(neighbours, basisVectors[0]) - 2 * v + GetNeighbourVelocity(neighbours, basisVectors[1]);
        c5 = (vis / Mathf.Pow(dy, 2)) * GetNeighbourVelocity(neighbours, basisVectors[2]) - 2 * v + GetNeighbourVelocity(neighbours, basisVectors[3]);
        c6 = (vis / Mathf.Pow(dz, 2)) * GetNeighbourVelocity(neighbours, basisVectors[4]) - 2 * v + GetNeighbourVelocity(neighbours, basisVectors[5]);
        return c1 + c2 + c3 + g + c4 + c5 + c6;
    }


    private float Divergence(Vector3Int pos)
    {
        float cx;
        float cy;
        float cz;
        Dictionary<Vector3Int, float> velocities;

        velocities = GetCellAllVelocities(gridArray[pos.x, pos.y, pos.z]);
        cx = (1 / dx) * (velocities[basisVectors[0]] - velocities[basisVectors[1]]);
        cy = (1 / dy) * (velocities[basisVectors[2]] - velocities[basisVectors[3]]);
        cz = (1 / dz) * (velocities[basisVectors[4]] - velocities[basisVectors[5]]);
        return -(cx + cy + cz);
    }

    private void VelocityUpdate(GridCell c)
    {

        float dp;
        Vector3Int b;
        Vector3Int pos = c.GetPos();

        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = GetCellAllNewVelocities(c); //Why is this new velocities when surface and solid use the current?
        float[] ds = new float[3] { dx, dy, dz };
        float d = Divergence(pos);
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
        cellsVelocities_X.SetNewVelocities(velocities, pos);
        cellsVelocities_Y.SetNewVelocities(velocities, pos);
        cellsVelocities_Z.SetNewVelocities(velocities, pos);

        c.SetNewPressure(c.GetPressure() + dp);
    }





// GOT THIS FAR BEFORE MY MIND WAS MELTING. 
// Apologies if code is wonky I tried my best
// TBH I shouldn't have attempted this alone as there are alot of design decisions that you might find weird apologies.
// Its 12am and i've got to try commiting this bad boy, so all the code below is commented out.
// All the original code is left in comments above the new functions I wrote (except Divergence whoops), so if it all fucks up then just comment my stuff out.
// Good luck :3


   private void SurfaceVelocityUpdate(GridCell c)
    {
        VelocityGrid[] velocityGrid = new VelocityGrid[3] {cellsVelocities_X, cellsVelocities_Y, cellsVelocities_Z};
        float p = 0f;
        Vector3Int pos = c.GetPos();
        Dictionary<Vector3Int, GridCell> neighbours = GetNeighbours(c);
        List<GridCell> inflowNeighbours = new List<GridCell>();
        List<GridCell> outflowNeighbours = new List<GridCell>();
        List<GridCell> openNeighbours = new List<GridCell>();
        List<Vector3Int> inflowDirections = new List<Vector3Int>();
        List<Vector3Int> outflowDirections = new List<Vector3Int>();
        List<Vector3Int> openDirections = new List<Vector3Int>();
        float totalInflow = 0;
        float totalOutflow = 0;
        int gridCount = -1;
        for (int i = 0; i < 6; i++)
        {
            if (i % 2 == 0)
            {
                gridCount++;
            }
            if (neighbours[basisVectors[i]].GetContents().Equals(Contents.Full) || neighbours[basisVectors[i]].GetContents().Equals(Contents.Surface))
            {
                if (basisVectors[i].x + basisVectors[i].y + basisVectors[i].z < 0)
                {
                    if (velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z] > 0)
                    {
                        inflowNeighbours.Add(neighbours[basisVectors[i]]);
                        inflowDirections.Add(basisVectors[i]);
                        totalInflow += velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z];
                    }
                    else
                    {
                        outflowNeighbours.Add(neighbours[basisVectors[i]]);
                        outflowDirections.Add(basisVectors[i]);
                        totalOutflow += velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z] * -1;
                    }
                }
                else
                {
                    pos = pos + basisVectors[i];
                    if (velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z] < 0)
                    {
                        inflowNeighbours.Add(neighbours[basisVectors[i]]);
                        inflowDirections.Add(basisVectors[i]);
                        totalInflow += velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z] * -1;
                    }
                    else
                    {
                        outflowNeighbours.Add(neighbours[basisVectors[i]]);
                        outflowDirections.Add(basisVectors[i]);
                        totalOutflow += velocityGrid[gridCount].GetVelocityGrid()[pos.x, pos.y, pos.z];
                    }
                    pos = pos - basisVectors[i];
                }
            }
            if (neighbours[basisVectors[i]].GetContents().Equals(Contents.Empty))
            {
                openNeighbours.Add(neighbours[basisVectors[i]]);
                openDirections.Add(basisVectors[i]);
            }
        }
        if (totalInflow - totalOutflow > 0)
        {
            foreach (Vector3Int d in openDirections)
            {   // May need to not use current velocity values
                if (d.x + d.y + d.z < 0)
                {
                    if (d.x != 0)
                    {
                        velocityGrid[0].SetNewVelocity(velocityGrid[0].GetVelocityGrid()[pos.x, pos.y, pos.z] - (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }else if(d.y != 0)
                    {
                        velocityGrid[1].SetNewVelocity(velocityGrid[1].GetVelocityGrid()[pos.x, pos.y, pos.z] - (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }
                    else if(d.z != 0)
                    {
                        velocityGrid[2].SetNewVelocity(velocityGrid[2].GetVelocityGrid()[pos.x, pos.y, pos.z] - (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }
                }
                else
                {
                    pos = pos + d;
                    if (d.x != 0)
                    {
                        velocityGrid[0].SetNewVelocity(velocityGrid[0].GetVelocityGrid()[pos.x, pos.y, pos.z] + (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }
                    else if (d.y != 0)
                    {
                        velocityGrid[1].SetNewVelocity(velocityGrid[1].GetVelocityGrid()[pos.x, pos.y, pos.z] + (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }
                    else if (d.z != 0)
                    {
                        velocityGrid[2].SetNewVelocity(velocityGrid[2].GetVelocityGrid()[pos.x, pos.y, pos.z] + (totalInflow - totalOutflow) / openNeighbours.Count, pos);
                    }
                    pos = pos - d;
                }
            }
        }
        c.SetNewPressure(p);
    }

    private void SolidVelocityUpdate(GridCell c)
    {
        float p = 0f;
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = GetCellAllVelocities(c);
        Dictionary<Vector3Int, GridCell> neighbours = GetNeighbours(c);
        Dictionary<Vector3Int, float> neighbourVelocities;
        GridCell neighbour;
        Vector3Int pos = c.GetPos();


        for (int i = 0; i < 6; i++)
        {
            if (neighbours.TryGetValue(basisVectors[i], out neighbour))
            {
                if (!neighbour.GetContents().Equals(Contents.Solid))
                {
                    neighbourVelocities = GetCellAllVelocities(neighbour);
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
            velocities = InstVelocities(new float[] { 0f, 0f, 0f, 0f, 0f, 0f});
        }
        cellsVelocities_X.SetNewVelocities(velocities, pos);
        cellsVelocities_Y.SetNewVelocities(velocities, pos);
        cellsVelocities_Z.SetNewVelocities(velocities, pos);
        c.SetNewPressure(p);
    }

    private void EmptyVelocityUpdate(GridCell c)
    {
        float p = 0f;
        Vector3Int pos = c.GetPos();
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities = GetCellAllVelocities(c);
        Dictionary<Vector3Int, GridCell> neighbours = GetNeighbours(c);
        GridCell neighbour;
        velocities = currentVelocities;
        for(int i = 0; i < 6; i++)
        {
            if (neighbours.TryGetValue(basisVectors[i], out neighbour))
            {
                velocities[basisVectors[i]] = GetCellAllVelocities(neighbour)[basisVectors[i] * -1];
            }
            else
            {

            }
        }
        cellsVelocities_X.SetNewVelocities(velocities, pos);
        cellsVelocities_Y.SetNewVelocities(velocities, pos);
        cellsVelocities_Z.SetNewVelocities(velocities, pos);
        c.SetNewPressure(p);
    }

    private void PressureIterations()
    {
        Dictionary<Vector3Int, float> velocities = new Dictionary<Vector3Int, float>();
        Dictionary<Vector3Int, float> currentVelocities;
        Vector3Int pos;
        Vector3Int b;
        GridCell c;
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
                    c = gridArray[x, y, z];
                    if (c.GetContents() == Contents.Full)
                    {
                        currentVelocities = GetCellAllVelocities(c);
                        pos = c.GetPos();
                        velocities.Clear();
                        for (int i = 0; i < 6; i++)
                        {
                            b = basisVectors[i];
                            velocities.Add(b, VelocityTilde(c, currentVelocities[b], b, GetNeighbours(c)));
                        }
                        cellsVelocities_X.SetNewVelocities(velocities, pos);
                        cellsVelocities_Y.SetNewVelocities(velocities, pos);
                        cellsVelocities_Z.SetNewVelocities(velocities, pos);
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
                        c = gridArray[x, y, z];
                        if (c.GetContents().Equals(Contents.Full))
                        {
                            VelocityUpdate(c);
                            if (c.GetDivergence() > eps)
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
                    c = gridArray[x, y, z];
                    if (c.GetContents() == Contents.Surface)
                    {
                        Debug.Log("Entered surface");
                        SurfaceVelocityUpdate(c);
                    }
                    else if (c.GetContents() == Contents.Solid)
                    {
                        SolidVelocityUpdate(c);
                    }
                    else if (c.GetContents() == Contents.Empty)
                    {
                        EmptyVelocityUpdate(c);
                    }
                    c.UpdateValues();
                }
            }
        }
        cellsVelocities_X.UpdateVelocities();
        cellsVelocities_Y.UpdateVelocities();
        cellsVelocities_Z.UpdateVelocities();
    }

    private VelocityDistance GetParticleVelocityComponents(Vector3 particlePosition, GridCell cell, Vector3 localCellPosition, Vector3 v1,Vector3Int v2,Vector3 v3,Vector3Int v4,Vector3 v5,Vector3Int v6)
    {
        Vector3[] velocityPositions = new Vector3[4];
        VelocityDistance componentVelocities = new VelocityDistance();
        Dictionary<Vector3Int, float> cellVelocities = GetCellAllVelocities(cell);
        componentVelocities.distance = new float[4];
        componentVelocities.velocity = new float[4];
        velocityPositions[0] = localCellPosition + (v1 * cellSize / 2);
        componentVelocities.velocity[0] = cellVelocities[v2];
        velocityPositions[1] = localCellPosition + (v3 * cellSize / 2);
        componentVelocities.velocity[1] = cellVelocities[v4];
        velocityPositions[2] = localCellPosition + (v5 * cellSize) + (v1 * cellSize / 2);
        componentVelocities.velocity[2] = GetCellAllVelocities(GetNeighbours(cell)[v6])[v2];
        velocityPositions[3] = localCellPosition + (v5 * cellSize) + (v3 * cellSize / 2);
        componentVelocities.velocity[3] = GetCellAllVelocities(GetNeighbours(cell)[v6])[v4];
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
        newPos = p.getPosition() + (v * dt);
        cellPos = water_grid.LocalToCell(newPos);
        if (cellPos[0] >= width-1 || cellPos[0] <= 0)
        {
            v[0] = 0f;
        }
        if (cellPos[1] >= height - 1 || cellPos[1] <= 0)
        {
            v[1] = 0f;
        }
        if (cellPos[2] >= depth - 1 || cellPos[2] <= 0)
        {
            v[2] = 0f;
        }
        newPos = p.getPosition() + (v * dt);
        cellPos = water_grid.LocalToCell(newPos);
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
        particleCount++;
        if (particleCount == 50)
        {
            particleCount = 0;
            CreateParticle(inflowLocation);
        }
        CellContents();
        Debug.Log("Update");
        PressureIterations();
        foreach (Particle p in particleList)
        {
            ParticleLocationUpdate(p);
        }
        //ParticleLocationUpdate(particle);
    }
}
