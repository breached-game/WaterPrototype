using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GridColumn : ScriptableObject
{
    private Vector2Int position;
    private int h;
    private int H;
    private int newh;
    private GameObject waterColumn;
    public void Setup(Vector2Int arg_position, int arg_H, GameObject column)
    {
        waterColumn = column;
        position = arg_position;
        H = arg_H;
        h = 0;
    }

    public Vector2Int GetPos()
    {
        return position;
    }

    public void UpdateValues()
    {
        h = newh;
    }

    public void SetNewh(int arg_h)
    {
        newh = arg_h;
    }

    public int Geth()
    {
        return h;
    }

    public void Seth(int arg_h)
    {
        h = arg_h;
    }

    public int GetH()
    {
        return H;
    }

    public void SetColumn()
    {
        waterColumn.transform.localScale = new Vector3(0.2f, 0.2f * h, 0.2f);
        waterColumn.transform.position = new Vector3(position.x * 0.2f, (0.2f * h) / 2, position.y * 0.2f);
        Debug.Log("Position: " + position);
        Debug.Log("Depth: " + h);
    }
}

