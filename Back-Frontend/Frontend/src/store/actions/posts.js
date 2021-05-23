import {
  RECEIVE_POSTS,
  REQUEST_POSTS,
  FAILURE_POSTS,
  RECEIVE_POST,
  REQUEST_POST,
  FAILURE_POST,
  UPDATE_POST_DATA,
  CLEAR_POST_DATA,
  CREATE_POST,
} from '../constants/posts';
import {
  postsList,
  getPostById as postById,
  createPost,
} from '../../services/posts';

//list of posts
const requestPosts = () => ({
  type: REQUEST_POSTS,
});

const receivePosts = (data) => ({
  type: RECEIVE_POSTS,
  data,
});

const failurePosts = (error) => ({
  type: FAILURE_POSTS,
  error,
});

// post by id

const requestPost = () => ({
  type: REQUEST_POST,
});

const receivePost = (data) => ({
  type: RECEIVE_POST,
  data,
});

const failurePost = (error) => ({
  type: FAILURE_POST,
  error,
});

export const updateData = (payload) => ({
  type: UPDATE_POST_DATA,
  payload,
});

export const clearPostData = () => ({
  type: CLEAR_POST_DATA,
});

export const newPost = (newpost) => async (dispatch) => {
  try {
    const data = await createPost(newpost);
    dispatch({
      type: CREATE_POST,
      payload: { userid: 1, ...data, id: new Date() },
    });
  } catch (error) {
    console.log(error);
  }
};

export const getPosts = () => async (dispatch) => {
  dispatch(requestPosts());
  try {
    const data = await postsList();
    dispatch(receivePosts(data));
  } catch (error) {
    dispatch(failurePosts(error));
  }
};

export const getPostById = (id) => async (dispatch) => {
  dispatch(requestPost());
  try {
    const data = await postById(id);
    dispatch(receivePost(data));
  } catch (error) {
    dispatch(failurePost(error));
  }
};
