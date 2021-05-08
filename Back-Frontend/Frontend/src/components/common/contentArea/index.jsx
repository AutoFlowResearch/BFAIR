import React, { Component } from "react";
import { Redirect, Route, Switch } from "react-router-dom";
// import { Route } from "../protectedRoutes";
import Dashboard from "../../pages/admin/dashboard";
import UserForm from "../../pages/admin/NewUser";
import Posts from "../../pages/admin/posts";
import NewPost from "../../pages/admin/posts/NewPost";
import NewRecipe from "../../pages/admin/recipe/NewRecipe";
import RecipeList from "../../pages/admin/recipe/RecipeList";
import SimpleForm from "../../pages/admin/SimpleForm";
import SampleSubmission from "../../pages/admin/userManagement"

class ContentArea extends Component {
  state = {
    role: "admin"
  };

  componentDidMount() {
    // console.log(this.props.location);
    // this.props.history.push(`/mainpage/${this.state.role}/dashboard`);
  }

  render() {
    return (
            
      <div className="container-fluid">
          <div className="row">
              <div className="col-lg-12 col-md-12 col-sm-12">
                  <div class="d-flex justify-content-between">
                      <b>Module</b>
                      <div class="main">
                          <div class="form-group has-search">
                              <span class="fa fa-search form-control-feedback"></span>
                              <input type="text" class="form-control" placeholder="Search for module name" />
                          </div>
                      </div>
                  </div>
              </div>
          </div>
          <div className="row">
              <div className="col-lg-6 col-md-6 col-sm-12">
                  <div className="card Mainview">
                      <div className="card-body">
                          <div className="row">
                              <div className="col-lg-8 col-md-8 col-sm-8">
                                  <div className="innerText">
                                      <b>DOE</b>
                                      <p>
                                          It has survived not only five centuries,
                                          but also the leap into electronic
                                      </p>
                                      <a href="bootstrap_flex.asp" className="readmore">Read More<i class="fas fa-arrow-right"></i></a>
                                  </div>
                              </div>
                              <div className="col-lg-4 col-md-4 col-sm-4">
                                  <div className="innerBox">
                                      <b>DOE</b>
                                  </div>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
              <div className="col-lg-6 col-md-6 col-sm-12">
                  <div className="card Mainview">
                      <div className="card-body">
                          <div className="row">
                              <div className="col-lg-8 col-md-8 col-sm-8">
                                  <div className="innerText">
                                      <b>Result Retrival</b>
                                      <p>
                                          It has survived not only five centuries,
                                          but also the leap into electronic
                                      </p>
                                      <a href="bootstrap_flex.asp" className="readmore">Read More<i class="fas fa-arrow-right"></i></a>
                                  </div>
                              </div>
                              <div className="col-lg-4 col-md-4 col-sm-4">
                                  <div className="innerBox">
                                      <b>Result Retrival</b>
                                  </div>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
          </div>
          <div className="row pt-3">
              <div className="col-lg-6 col-md-6 col-sm-12">
                  <div className="card Mainview">
                      <div className="card-body">
                          <div className="row">
                              <div className="col-lg-8 col-md-8 col-sm-8">
                                  <div className="innerText">
                                      <b>Module 3</b>
                                      <p>
                                          It has survived not only five centuries,
                                          but also the leap into electronic
                                      </p>
                                      <a href="bootstrap_flex.asp" className="readmore">Read More<i class="fas fa-arrow-right"></i></a>
                                  </div>
                              </div>
                              <div className="col-lg-4 col-md-4 col-sm-4">
                                  <div className="innerBox">
                                  <b>Module 3</b>
                                  </div>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
              <div className="col-lg-6 col-md-6 col-sm-12">
                  <div className="card Mainview">
                      <div className="card-body">
                          <div className="row">
                              <div className="col-lg-8col-md-8 col-sm-8">
                                  <div className="innerText">
                                      <b>Module 4</b>
                                      <p>
                                          It has survived not only five centuries,
                                          but also the leap into electronic
                                      </p>
                                      <a href="bootstrap_flex.asp" className="readmore">Read More<i class="fas fa-arrow-right"></i></a>
                                  </div>
                              </div>
                              <div className="col-lg-4 col-md-4 col-sm-4">
                                  <div className="innerBox">
                                  <b>Module 4</b>
                                  </div>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
          </div>
          <div className="row mt-3">
              <div className="col-lg-12 col-md-12 col-sm-12">
                  <h5>Project Source Code</h5>
                  <div className="card bottomCard">
                      <div className="card-body">
                          <div className="row">
                              <div className="col-lg-2 col-md-2 col-sm-2">
                                  <div className="innerBox">Lawrence</div>
                              </div>
                              <div className="col-lg-10 col-md-10 col-sm-10">
                                  <div className="innerText">
                                      <b>Source code</b>
                                      <p>
                                          It has survived not only five centuries,
                                          but also the leap into electronic
                                      </p>
                                      <a href="bootstrap_flex.asp" className="readmore">Read More<i class="fas fa-arrow-right"></i></a>
                                  </div>
                              </div>
                          </div>
                      </div>
                  </div>
              </div>
          </div>
      </div>

  )
  }
}

export default ContentArea;
