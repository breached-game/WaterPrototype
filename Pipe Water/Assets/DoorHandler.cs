using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DoorHandler : MonoBehaviour
{
    public GameObject door;
    public GameObject mainWall;
    public GameObject wall1;
    public GameObject wall2;
    public GameObject openDoorButton;
    public GameObject closeDoorButton;
    public GameObject waterGrid;

    public void OpenDoor()
    {
        door.SetActive(false);
        waterGrid.GetComponent<WaterGrid>().RemoveObject(door);
        openDoorButton.SetActive(false);
        closeDoorButton.SetActive(true);
    }

    public void CloseDoor()
    {
        mainWall.SetActive(true);
        wall1.SetActive(false);
        wall2.SetActive(false);
        openDoorButton.SetActive(true);
        closeDoorButton.SetActive(false);
    }
}
