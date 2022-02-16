using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GridColumn : ScriptableObject
{
    private Vector2Int position;
    private float h;
    private float H;
    private float newh;
    private GameObject waterColumn;
    private Dictionary<Vector2Int, float> newOutflows;
    private Dictionary<Vector2Int, float> outflows;
    public void Setup(Vector2Int arg_position, float arg_H, GameObject column)
    {
        waterColumn = column;
        position = arg_position;
        H = arg_H;
        h = 0f;
        outflows = new Dictionary<Vector2Int, float>();
        outflows.Add(Vector2Int.right, 0f);
        outflows.Add(Vector2Int.left, 0f);
        outflows.Add(Vector2Int.up, 0f);
        outflows.Add(Vector2Int.down, 0f);
        newOutflows = new Dictionary<Vector2Int, float>();
        newOutflows.Add(Vector2Int.right, 0f);
        newOutflows.Add(Vector2Int.left, 0f);
        newOutflows.Add(Vector2Int.up, 0f);
        newOutflows.Add(Vector2Int.down, 0f);
    }

    public Vector2Int GetPos()
    {
        return position;
    }

    public void UpdateValues()
    {
        h = newh;
        outflows = newOutflows;
    }

    public void SetNewh(float arg_h)
    {
        newh = arg_h;
    }

    public float Geth()
    {
        return h;
    }

    public void Seth(float arg_h)
    {
        h = arg_h;
    }

    public float GetH()
    {
        return H;
    }

    public void SetColumn()
    {
        if (h < 0)
        {
            Debug.Log(position);
            Debug.Log("Depth: " + h);
        }
        waterColumn.transform.localScale = new Vector3(0.2f, 0.2f * h, 0.2f);
        waterColumn.transform.position = new Vector3(position.x * 0.2f, (0.2f * (h + H)) / 2, position.y * 0.2f);
        //Debug.Log("Position: " + position);
        //Debug.Log("Depth: " + h);
    }

    public Dictionary<Vector2Int, float> GetNewOutflows()
    {
        return newOutflows;
    }

    public Dictionary<Vector2Int, float> GetOutflows()
    {
        return outflows;
    }

    public void SetNewOutflows(Dictionary<Vector2Int, float> fs)
    {
        newOutflows[Vector2Int.right] = fs[Vector2Int.right];
        newOutflows[Vector2Int.left] = fs[Vector2Int.left];
        newOutflows[Vector2Int.up] = fs[Vector2Int.up];
        newOutflows[Vector2Int.down] = fs[Vector2Int.down];
    }

    public void SetOutflows(Dictionary<Vector2Int, float> fs)
    {
        outflows = fs;
    }
}
