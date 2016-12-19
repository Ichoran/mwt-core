resolvers += "Sonatype OSS Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots"

lazy val root = (project in file(".")).
  settings(
    scalaVersion := "2.12.1",
    name := "mwt-scala-wcon",
    version := "0.1.0",
    scalacOptions ++= Seq("-unchecked", "-feature", "-deprecation"),
    libraryDependencies += "com.github.ichoran" %% "kse" % "0.6-SNAPSHOT"
  )
