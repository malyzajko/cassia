import sbt._
import Keys._

object BuildSettings {
  val buildSettings = Defaults.defaultSettings ++ Seq(
    organization := "ch.epfl.lara",
    version := "1.0.0",
    scalaVersion := "2.10.2",
    resolvers += Resolver.sonatypeRepo("snapshots"),
    resolvers += Resolver.sonatypeRepo("releases"),
    scalacOptions ++= Seq()
  )
}

object MyBuild extends Build {
  import BuildSettings._

  lazy val root: Project = Project(
    "root",
    file("."),
    settings = buildSettings ++ Seq(
      run <<= run in Compile in core)
  ) aggregate(macros, core)

  lazy val macros: Project = Project(
    "macros",
    file("macros"),
    settings = buildSettings ++ Seq(
      libraryDependencies <+= (scalaVersion)("org.scala-lang" % "scala-reflect" % _))
  )

  lazy val core: Project = Project(
    "applications",
    file("applications"),
    settings = buildSettings
  ) dependsOn(macros)
}
