## Apriltags detection written in C#, based on Pupil Labs C library.

A translation of the pupil labs library for apriltags, very slight variation in results. 
Unlimited robotics allows for apriltags detection inside the simulation of our robot, Gary, for a number of tasks such as inspecting inventories or marking world positions for the robot to place objects at. We are happy to share the code with anyone else who wishes to use this great tool in his project. It is very important for us to emphasize that this would not have been possible without Pupil Labs amazing product and that most of the credit lies with them.

https://github.com/pupil-labs/apriltags

The library is built for unity, to get rid of unity dependencies(very small amounts of debug code) you have to edit two scripts, remove anything between the comments "// unity code start" and "// unity code end" if you desire to use a regular C# solution. 

Scripts with unity dependencies: ApriltagLog, ApriltagPrint, ApriltagTest, ApriltagImage
