from pydantic import BaseModel, Field


# Target box information for target
class TargetBox(BaseModel):
    center_x: int
    center_y: int
    center_z: int
    size_x: int
    size_y: int
    size_z: int


# Request and response data models
class DockingRequest(BaseModel):
    bucket_name: str = Field(..., "S3 bucket name")
    selfies_object_name: str = Field(..., "S3 object name to smi file")
    target_object_name: str = Field(..., "S3 object name for target pdb file")
    target_box: TargetBox = Field(..., "Object containing bounding box dimensions")


class DockingResponse(BaseModel):
    pass
    